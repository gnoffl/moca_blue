# moca_blue
[2023-06-07]

MOCA BLUE

MOtif
  Characterization
&  Annotation
   from DEEP LEARNING feature enrichment

Welcome the the moca_blue suite!
from Simon M. Zumkeller
RStudio
2022.07.2 Build 576

This is a tool-box for the analyses of DNA motifs
that have been derived from deep-learning model features extraction.
moca_blue is currently in development.



---------- Workflow example --------
mkdir moca_blue
cd mocablue

mkdir 0MOTIFS mo_nom mo_range mo_proj mo_clu ref_seq ...

./mo_nom/rdf5_get_cwms_per_pattern.v1.0.R

./mo_clu/Mo_cluster_v2.0.R

./mo_range/rdf5_get_seql_patternV2.1.R
./mo_range/meta_motif_ranges_characteristics_TSS-TTS.1.4.R

#Map motifs.jaspar to reference genome using BLAMM (https://github.com/biointec/blamm)

./blamm-master/build/blamm_meV1.0.sh
cp ./blamm-master/build/outPROJECT ./moca_blue/mo_proj/outPROJECT

./mo_proj/occ_filter_v1.1R
./mo_proj/mo_feat-filt.v2.3.R
./mo_proj/mo_feature_tester.v1.0.R
./mo_proj/mo_predictabilityV1.0.R






---------- Overview --------
recommended directory structure
Please find more detailed descriptions
of the directories and the code in the files itself.

This is a pipeline of consecutive operations that can be and will be availabe here.

INPUT DIRECTORY                /0MOTIFS                              /ref_seq
                      - HDF5.file [feature extraction files]       - fastas
                      _____|_________________                      - gffs
                      |                      |                     - meta-data
START DIRECTORY   /mo_nom                   /mo_range                |
output          - get motif patterns      - get motif meta-data      |
                - motif annotation           |                       |
                - motif modification                                 |
                      |______________________________________________|
                      |                                |
                  /mo_clu                             MAPPING to reference (external)
                - analyze motifs             |        use e.g. "blamm
                - compare/cluster            |        (https://github.com/biointec/blamm)
                                             |        cp occurences.txt [results] /mo_proj
                                             |_________|
                                                  |
                                                 /mo_proj
                                                - filter for meaningful matches
                                                - interpret model predictions
                                                - gene annotation
                                                - module generation
                                                - get bed files for visualisation in JBROWSE



mo_nom  --------------------------

Extract motifs from MoDisco hdf5 files and assign nomenclature.
Currently, there are three versions of the same script that can be used for the extraction of a given format of weight matrix.

rdf5_get_xxx_per_pattern.v1.0R.R

PFM - positional frequence matrix
PWM - positional weight matrix (best for clustering/comparison)
CWM - contribution weight matrix (best for mapping)

mo_range ------------------------

Motifs/ EPMs are not distributed at random in a genome.
To optimize the search for motifs/EPMs in a genome or gene-space, these tools
extract the positionally preferred ranges for each motif/EPM in a hdf5 file.

rdf5_get_seql_per_patternV2.R - Extract a list of seqlets and their positions from the hdf5 file

meta_motif_ranges_characteristics_TSS-TTS.1.1.R - Producee a table from the rdf5_get_seql_per_patternV2.R output
  that provides the gene-space statistics for each motif/seqlet in reference to transcription start and stop sites (TSS, TTS)


mo_clu --------------------------

Analyse and Edit motif-files stored in jaspar-format here. Results should be stored in the "out" directory.

mo_cluster_v2.0R - generates dendrograms/trees based on distancy-matrix for different models.

mo_old ---------------------------

Old and outdated scripts used for the moca_blue suite are stored here.


ref_seq -------------------------

Store genome data like fasta, gff and many more here for INPUT. 

mo_proj -------------------------
After searching for motif matches within a reference sequence, these results can be tested and further characterized using the occ_filter, mo_feat-filter, mo_feature_tester and mo_predictability scripts.

