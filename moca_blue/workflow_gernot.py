import os
import argparse
import json
from pydantic import BaseModel, field_validator
from typing import Dict, List, Optional, Tuple
import subprocess
import shutil


class Config(BaseModel):
    base_name: str
    hdf5_path: str
    fasta_path: str

    @field_validator("base_name", mode='after')
    def validate_base_name(cls, v: str) -> str:
        if not v:
            raise ValueError("Base name cannot be empty.")
        return v

    @field_validator("hdf5_path", "fasta_path", mode='after')
    def validate_paths(cls, v: str) -> str:
        if not os.path.exists(v):
            raise ValueError(f"File {v} does not exist.")
        return v


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="mocaBlueAnnotationAndFiltering", description="Performs the full process from extracting motifs to" \
    "mapping and filtering them based on user-defined criteria.")
    parser.add_argument("--input", "-i", type=str, required=True, help="Path to the input/config file")
    args = parser.parse_args()
    return args


def load_configs(input_path) -> Tuple[List[Config], List[Tuple[str, str]]]:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file {input_path} does not exist.")

    with open(input_path, 'r') as file:
        configs = json.load(file)
    
    working_configs, errors = [], []
    for config in configs:
        try:
            config_model = Config(**config)
            working_configs.append(config_model)
        except ValueError as e:
            errors.append((json.dumps(config, indent=2), str(e)))

    return working_configs, errors


def run_command(*args: str, output: bool = False) -> Optional[str]:
    print(f"Running command: {' '.join(args)}")
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        raise RuntimeError(f"Command failed with error message: {result.stderr}")
    if output:
        return result.stdout


def extract_motifs(base_name: str, hdf5_path: str, output_base_dir: str) -> str:
    print(f"Extracting motifs for {base_name} from {hdf5_path}...")
    output_folder = os.path.join(output_base_dir, "motifs")
    os.mkdir(output_folder)
    run_command("Rscript", "moca_blue/mo_nom/get_rdf5_cwms_per_pattern_v1.1.R", ".", base_name, hdf5_path, output_folder)
    print(f"Motifs extracted to {output_folder}.")
    return output_folder


def map_motifs(motifs_path: str, fasta_path: str, output_base_dir: str, basename: str) -> str:
    output_folder = os.path.join(output_base_dir, "mapped_motifs")
    motif_files = [os.path.join(motifs_path, f) for f in os.listdir(motifs_path) if f.endswith('.jaspar')]
    if len(motif_files) != 1:
        raise ValueError(f"Expected exactly one file in {motifs_path}, found {len(motif_files)}.")
    full_motif_path = motif_files[0]
    os.mkdir(output_folder)
    print(f"writing multi fasta file to {output_folder}")
    mf_path = os.path.join(output_folder, "sequences.mf")
    with open(mf_path, 'w') as multi_fasta:
        multi_fasta.write(f"{basename}\t{fasta_path}\n")
    print(f"Mapping motifs from {motifs_path} to {fasta_path}...")
    run_command("./blamm_meV1.0.sh", output_folder, full_motif_path, mf_path)
    print(f"mappings saved to {output_folder}.")
    return output_folder


def calculate_ranges(hdf5_path: str, output_base_dir: str, base_name: str) -> str:
    output_folder = os.path.join(output_base_dir, "ranges")
    os.mkdir(output_folder)
    print(f"Extracting seqlets from {hdf5_path}...")
    run_command("Rscript", "moca_blue/mo_ran/rdf5_get_seql_per_patternV2.1.R", ".", base_name, hdf5_path, output_folder)
    print(f"seqlets extracted to {output_folder}.")
    seqlet_files = [os.path.join(output_folder, f) for f in os.listdir(output_folder) if f.endswith('.txt')]
    if len(seqlet_files) != 1:
        raise ValueError(f"Expected exactly one seqlet file in {output_folder}, found {len(seqlet_files)}.")
    seqlet_path = seqlet_files[0]
    print(f"Calculating ranges for seqlets in {seqlet_path}...")
    run_command("Rscript", "moca_blue/mo_ran/epm_occurence_ranges_TSS-TTS.1.6.R", ".", seqlet_path, hdf5_path, base_name, output_folder)
    print(f"Ranges calculated and saved to {output_folder}.")
    return output_folder


def create_gff(fasta_path: str) -> str:
    gff_path = os.path.splitext(fasta_path)[0] + ".gff3"
    if os.path.exists(gff_path):
        print(f"found annotation file {gff_path}")
        return gff_path
    print(f"Creating GFF file from {fasta_path}...")
    output = run_command("python", "moca_blue/ref_seq/create_annotation.py", "-f", fasta_path, output=True)
    if output is None:
        raise RuntimeError("really shouldnt come to this. run_command somehow didnt return output :(")
    last_line = output.strip().split('\n')[-1]
    path = last_line.split()[-1]
    print(f"GFF file created at {path}")
    return path


def filter_mapping(mapping_path: str, ranges_path: str, output_base_dir: str, base_name: str, gff_path: str) -> None:
    print(f"Filtering mapping results from {mapping_path}...")
    output_folder = os.path.join(output_base_dir, "filtered")
    os.mkdir(output_folder)
    simple_filtered_mapping_path = os.path.join(output_folder, f"{base_name}_filtered.tsv")
    mapping_file_path = os.path.join(mapping_path, "occurrences.txt")
    run_command("Rscript", "moca_blue/mo_proj/occ_filter_v1.1.R", ".", mapping_file_path, simple_filtered_mapping_path)
    print(f"Initial filtered mapping results saved to {simple_filtered_mapping_path}")
    ranges_files = [os.path.join(ranges_path, f) for f in os.listdir(ranges_path)]
    tss_paths = [f for f in ranges_files if f.endswith("TSS_motif_ranges_q1q9.csv")]
    tts_paths = [f for f in ranges_files if f.endswith("TTS_motif_ranges_q1q9.csv")]
    if len(tss_paths) != 1 or len(tts_paths) != 1:
        raise ValueError(f"Expected exactly one TSS and TTS file in {ranges_path}, found {len(tss_paths)} TSS and {len(tts_paths)} TTS files.")
    tss_path = tss_paths[0]
    tts_path = tts_paths[0]
    print(f"Applying range filters from {ranges_path}...")
    run_command("Rscript", "moca_blue/mo_proj/mo_feat-filter.v3.4.R", ".", base_name, output_folder, tss_path, tts_path, gff_path, simple_filtered_mapping_path)
    print(f"Filtered mapping results saved to {output_folder}")


def run_workflow(config: Config):
    print(f"{'':_>50}")
    print(f"Running workflow for base_name: {config.base_name}, hdf5_path: {config.hdf5_path}, fasta_path: {config.fasta_path}")
    output_base_dir = f"moca_blue/out/{config.base_name}"
    if os.path.exists(output_base_dir):
        raise FileExistsError(f"Output directory {output_base_dir} already exists. Please remove it or choose a different base name.")
    os.makedirs(output_base_dir, exist_ok=True)
            
    motifs_path = extract_motifs(config.base_name, config.hdf5_path, output_base_dir=output_base_dir)
    ranges_path = calculate_ranges(hdf5_path=config.hdf5_path, output_base_dir=output_base_dir, base_name=config.base_name)
    mapping_path = map_motifs(motifs_path=motifs_path, fasta_path=config.fasta_path, output_base_dir=output_base_dir, basename=config.base_name)
    gff_path = create_gff(config.fasta_path)
    filter_mapping(mapping_path, ranges_path, output_base_dir=output_base_dir, base_name=config.base_name, gff_path=gff_path)


def print_errors(errors):
    print(f"\n{'RUN RESULTS':_^50}\n")
    if errors:
        print("Errors in configuration:")
        for error in errors:
            print(f"{'':_>50}\n")
            print(f"Config:\n{error[0]}\n\nError:\n{error[1]}")
            print(f"{'':_>50}\n")
    else:
        print(f"{'':_>50}\n")
        print("all configs ran successfully!")


def main():
    args = parse_args()
    configs, errors = load_configs(args.input)
    
    for config in configs:
        try:
            run_workflow(config)
        except Exception as e:
            errors.append((json.dumps(config.model_dump(), indent=2), str(e)))

    print_errors(errors)


if __name__ == "__main__":
    main()