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

Please find more detailed descriptions
of the directories and their role within them, respectively.

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
