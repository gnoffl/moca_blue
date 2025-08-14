import os
import argparse
import json
import random
import pandas as pd
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


def extract_motifs(base_name: str, hdf5_path: str, output_base_dir: str) -> Tuple[str, bool]:
    output_folder = os.path.join(output_base_dir, "motifs")
    expected_output = os.path.join(output_folder, f"{base_name}_cwm-motifs.jaspar")
    
    if os.path.exists(expected_output):
        print(f"Motifs already exist at {output_folder}, skipping extraction.")
        return output_folder, True
    
    print(f"Extracting motifs for {base_name} from {hdf5_path}...")
    os.makedirs(output_folder, exist_ok=True)
    run_command("Rscript", "moca_blue/mo_nom/get_rdf5_cwms_per_pattern_v1.1.R", ".", base_name, hdf5_path, output_folder)
    print(f"Motifs extracted to {output_folder}.")
    return output_folder, False


def map_motifs(motifs_path: str, fasta_path: str, output_base_dir: str, basename: str) -> Tuple[str, bool]:
    output_folder = os.path.join(output_base_dir, "mapped_motifs")
    expected_output = os.path.join(output_folder, "occurrences.txt")
    
    if os.path.exists(expected_output):
        print(f"Mappings already exist at {output_folder}, skipping mapping.")
        return output_folder, True
    
    motif_files = [os.path.join(motifs_path, f) for f in os.listdir(motifs_path) if f.endswith('.jaspar')]
    if len(motif_files) != 1:
        raise ValueError(f"Expected exactly one file in {motifs_path}, found {len(motif_files)}.")
    full_motif_path = motif_files[0]
    os.makedirs(output_folder, exist_ok=True)
    print(f"writing multi fasta file to {output_folder}")
    mf_path = os.path.join(output_folder, "sequences.mf")
    with open(mf_path, 'w') as multi_fasta:
        multi_fasta.write(f"{basename}\t{fasta_path}\n")
    print(f"Mapping motifs from {motifs_path} to {fasta_path}...")
    run_command("./blamm_meV1.0.sh", output_folder, full_motif_path, mf_path)
    print(f"mappings saved to {output_folder}.")
    return output_folder, False


def calculate_ranges(hdf5_path: str, output_base_dir: str, base_name: str) -> Tuple[str, bool]:
    output_folder = os.path.join(output_base_dir, "ranges")
    expected_tss = os.path.join(output_folder, f"{base_name}-TSS_motif_ranges_q1q9.csv")
    expected_tts = os.path.join(output_folder, f"{base_name}-TTS_motif_ranges_q1q9.csv")
    
    if os.path.exists(expected_tss) and os.path.exists(expected_tts):
        print(f"Ranges already exist at {output_folder}, skipping range calculation.")
        return output_folder, True
    
    os.makedirs(output_folder, exist_ok=True)
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
    return output_folder, False


def create_gff(fasta_path: str) -> Tuple[str, bool]:
    gff_path = os.path.splitext(fasta_path)[0] + ".gff3"
    if os.path.exists(gff_path):
        print(f"found annotation file {gff_path}")
        return gff_path, True
    print(f"Creating GFF file from {fasta_path}...")
    output = run_command("python", "moca_blue/ref_seq/create_annotation.py", "-f", fasta_path, output=True)
    if output is None:
        raise RuntimeError("really shouldnt come to this. run_command somehow didnt return output :(")
    last_line = output.strip().split('\n')[-1]
    path = last_line.split()[-1]
    print(f"GFF file created at {path}")
    return path, False


def filter_mapping(mapping_path: str, ranges_path: str, output_base_dir: str, base_name: str, gff_path: str) -> Tuple[str, bool]:
    output_folder = os.path.join(output_base_dir, "filtered")
    expected_outputs = [
        os.path.join(output_folder, f"{base_name}_gene_none-q1q9.csv"),
        os.path.join(output_folder, f"{base_name}_gene_none-mima.csv"),
        os.path.join(output_folder, f"{base_name}_filtered.tsv")
    ]
    
    if all(os.path.exists(f) for f in expected_outputs):
        print(f"Filtered results already exist at {output_folder}, skipping filtering.")
        return output_folder, True

    print(f"Filtering mapping results from {mapping_path}...")
    os.makedirs(output_folder, exist_ok=True)
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
    return output_folder, False

def compare_motifs_jaspar(motifs_path: str, output_base_dir: str) -> Tuple[str, bool]:
    output_folder = os.path.join(output_base_dir, "jaspar_comparison")
    motif_files = [os.path.join(motifs_path, f) for f in os.listdir(motifs_path) if f.endswith('.jaspar')]
    if len(motif_files) != 1:
        raise ValueError(f"Expected exactly one motif file in {motifs_path}, found {len(motif_files)}.")
    
    motif_basename = os.path.splitext(os.path.basename(motif_files[0]))[0]
    expected_output = os.path.join(output_folder, f"{motif_basename}_comparison_JASPAR2020.csv")
    
    if os.path.exists(expected_output):
        print(f"JASPAR comparison already exists at {output_folder}, skipping comparison.")
        return expected_output, True
    
    print(f"Comparing motifs against JASPAR2020 database...")
    os.makedirs(output_folder, exist_ok=True)
    full_motif_path = motif_files[0]
    
    run_command("Rscript", "moca_blue/mo_nom/mo_compare_JASPAR2020_v1.0.R", ".", full_motif_path, output_folder)
    print(f"JASPAR comparison results saved to {output_folder}")
    return expected_output, False


def generate_motif_logos(motifs_path: str, output_base_dir: str, base_name: str) -> Tuple[str, bool]:
    """
    Generate sequence logos for motifs using the mo_cluster_v2.7.R script.
    Note: This version only generates PNG logos, clustering functionality is commented out.
    
    Args:
        motifs_path: Path to directory containing JASPAR motif files
        output_base_dir: Base output directory for the workflow
        base_name: Base name for the analysis
    
    Returns:
        Tuple of (output_folder_path, was_skipped)
    """
    output_folder = os.path.join(output_base_dir, "motif_logos")
    motif_files = [os.path.join(motifs_path, f) for f in os.listdir(motifs_path) if f.endswith('.jaspar')]
    if len(motif_files) != 1:
        raise ValueError(f"Expected exactly one motif file in {motifs_path}, found {len(motif_files)}.")
    
    # Check if any PNG files exist in the output folder (indicates logos were generated)
    if os.path.exists(output_folder) and any(f.endswith('.png') for f in os.listdir(output_folder)):
        print(f"Motif logos already exist at {output_folder}, skipping logo generation.")
        return output_folder, True
    
    print(f"Generating motif logos from {motifs_path}...")
    os.makedirs(output_folder, exist_ok=True)
    full_motif_path = motif_files[0]
    
    # Run motif logo generation script
    run_command("Rscript", "moca_blue/mo_clu/mo_cluster_v2.7.R", ".", full_motif_path, output_folder)
    print(f"Motif logos saved to {output_folder}")
    return output_folder, False


def update_transcription_factor_names(filtered_path: str, epm_tf_mapping_path: str) -> Tuple[str, bool]:
    """
    Update the mapping file to use transcription factor names instead of EPM names based on JASPAR comparison results.
    
    Args:
        mapping_path: Path to the directory containing the mapping results
        jaspar_path: Path to the JASPAR comparison results
    """
    skipped = False
    for ending in ["gene_none-q1q9.csv", "gene_none-mima.csv"]:
        tf_mapping_df = pd.read_csv(epm_tf_mapping_path, sep="\t")
        current_filtered = [os.path.join(filtered_path, file) for file in os.listdir(filtered_path) if file.endswith(ending)]
        if len(current_filtered) != 1:
            raise ValueError(f"Expected exactly one filtered file ending with {ending} in {filtered_path}, found {len(current_filtered)}.")
        filtered_df = pd.read_csv(current_filtered[0])
        if "tf_name" in filtered_df.columns:
            print(f"Transcription factor names already exist in {current_filtered[0]}, skipping update.")
            skipped = True
            continue
        # cleave end of epms
        tf_mapping_df['subj'] = tf_mapping_df['subj'].str.split('_').str[:-1].str.join('_').str[:-1]
        tf_mapping_df = tf_mapping_df.drop_duplicates(subset=['subj'])
        tf_mapping_df = tf_mapping_df[['subj', 'targ']].rename(columns={'subj': 'epm', 'targ': 'tf_name'})
        filtered_df = filtered_df.merge(tf_mapping_df, left_on='epm', right_on='epm', how='left')
        # change order of columns such that tf_name is the second column
        cols = filtered_df.columns.tolist()
        cols.remove('tf_name')
        cols.insert(1, 'tf_name')
        filtered_df = filtered_df[cols]
        filtered_df.to_csv(current_filtered[0], index=False)
        print(f"Updated transcription factor names in {current_filtered[0]}")
    return filtered_path, skipped


def run_workflow(config: Config):
    print(f"{'':_>50}")
    print(f"Running workflow for base_name: {config.base_name}, hdf5_path: {config.hdf5_path}, fasta_path: {config.fasta_path}")
    output_base_dir = f"moca_blue/out/{config.base_name}"
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Track which steps were skipped/executed
    step_status = {}
    
    motifs_path, motifs_skipped = extract_motifs(config.base_name, config.hdf5_path, output_base_dir=output_base_dir)
    step_status["extract_motifs"] = "skipped" if motifs_skipped else "executed"
    
    epm_tf_mapping_path, jaspar_skipped = compare_motifs_jaspar(motifs_path=motifs_path, output_base_dir=output_base_dir)
    step_status["compare_motifs_jaspar"] = "skipped" if jaspar_skipped else "executed"
    
    logos_path, logos_skipped = generate_motif_logos(motifs_path=motifs_path, output_base_dir=output_base_dir, base_name=config.base_name)
    step_status["generate_motif_logos"] = "skipped" if logos_skipped else "executed"
    
    ranges_path, ranges_skipped = calculate_ranges(hdf5_path=config.hdf5_path, output_base_dir=output_base_dir, base_name=config.base_name)
    step_status["calculate_ranges"] = "skipped" if ranges_skipped else "executed"
    
    mapping_path, mapping_skipped = map_motifs(motifs_path=motifs_path, fasta_path=config.fasta_path, output_base_dir=output_base_dir, basename=config.base_name)
    step_status["map_motifs"] = "skipped" if mapping_skipped else "executed"
    
    gff_path, gff_skipped = create_gff(config.fasta_path)

    filtered_path, filtering_skipped = filter_mapping(mapping_path, ranges_path, output_base_dir=output_base_dir, base_name=config.base_name, gff_path=gff_path)
    step_status["filter_mapping"] = "skipped" if filtering_skipped else "executed"

    # use results from motif comparison to use transcription factor names instead of epm names
    skipped_folders, skipped_tf_names = update_transcription_factor_names(filtered_path=filtered_path, epm_tf_mapping_path=epm_tf_mapping_path)
    step_status["update_transcription_factor_names"] = "skipped" if skipped_tf_names else "executed"

    # Generate summary
    executed_steps = [step for step, status in step_status.items() if status == "executed"]
    skipped_steps = [step for step, status in step_status.items() if status == "skipped"]
    
    if skipped_steps:
        error_message = "\n".join([f" - {step}" for step in skipped_steps])
        raise Warning(f"skipped the following steps, since the expected output was already present in {output_base_dir}:\n{error_message}")

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