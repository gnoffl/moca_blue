# moca_blue
[2024-03-13]  
Welcome the the moca_blue suite!  
from Simon M. Zumkeller  
RStudio  
2022.07.2 Build 576

This is a tool-box for the analyses of DNA motifs  
that have been derived from deep-learning model features extraction software 'deepCRE'.  
https://github.com/NAMlab/DeepCRE  
moca_blue is currently in development.  

################## ---------- Workflow example -------- ###################################  
#moca_blue has a defined directory structure that needs to be established for proper usage#  
mkdir moca_blue  
cd mocablue  
#Generate folder structure  
mkdir 0MOTIFS mo_nom mo_range mo_proj mo_clu ref_seq ..  

#Results will be stored in the out directories   
mkdir 0MOTIFS mo_nom/out mo_range/out mo_proj/out ..  

# Assign nomenclature and extract EPMs  into jaspar file format     #######################  
mo_nom
Extract motifs from TF-MoDisco hdf5 files, assigns nomenclature and produces weblogos.  
Currently, there are three versions of the same script that can be used for the extraction  
of a given format of weight matrix.  
PFM - positional frequence matrix  
PWM - positional weight matrix (best for clustering/comparison)  
CWM - contribution weight matrix (best for mapping)  

./mo_nom/get_rdf5_cwms_per_pattern_v1.1.R  

#Compare EPMs to JASPAR2020 database  
./mo_nom/mo_compare_JASPAR2020_v1.0.R  

# Compare EPMs to themselves, other files in jaspar format or JASPAR2020 database        ####  
mo_clu   

Analyse and Edit motif-files stored in jaspar-format here. Results should be stored in the "out" directory.  
Generates dendrograms/trees based on similarity-matrix for EPMs and Visual.  

./mo_clu/mo_cluster_v2.7.R  
./mo_clu/mo_compare_JASPAR2020_v1.0.R  

#Visualize EPM clustering results  
./mo_clu/mo_tree_viz.SZ.v1.0.R  

# Extract saliency maps and importace scores      ###########################################  
./mo_imp/rdf5_get_epm_contrib_scores.v1.1.R  
./mo_imp/mo_imp_scores.v1.1.R  

#Visualize saliency maps and importace scores (in development)  
./mo_imp/mo_imp_depth_v0.7.R  

# Extract seqlet occurring ranges, positional preferences of EPMs        ####################  
mo_range ------------------------  

Motifs/ EPMs are not distributed at random in a genome.  
To optimize the search for motifs/EPMs in a genome or gene-space, these tools  
extract the positionally preferred ranges for each motif/EPM in a hdf5 file.  

rdf5_get_seql_per_patternV2.R - Extract a list of seqlets and their positions from the hdf5 file  
meta_motif_ranges_characteristics_TSS-TTS.1.1.R - Producee a table from the rdf5_get_seql_per_patternV2.R output  
  that provides the gene-space statistics for each motif/seqlet in reference to transcription start and stop sites (TSS, TTS)  
  
./mo_ran/rdf5_get_seql_per_patternV2.1.R  
./mo_ran/epm_occurence_ranges_TSS-TTS.1.6.R  

mo_nom/rdf5_get_cwms_per_pattern.v1.0.R  
./mo_clu/Mo_cluster_v2.0.R  

./mo_range/rdf5_get_seql_patternV2.1.R  
./mo_range/meta_motif_ranges_characteristics_TSS-TTS.1.4.R  

# Map motifs.jaspar to reference genome using BLAMM (https://github.com/biointec/blamm) #####  
#Follow BLAMM installation guide   
mv ./blamm_meV1.0.sh ../blamm-master/build/blamm_meV1.0.sh  
#Edit sequences.mf file to specify target file for EPM search; e.g. ./ref_seq/file.fas  
./blamm-master/build/blamm_meV1.0.sh  
cp ./blamm-master/build/outPROJECT ./moca_blue/mo_proj/outPROJECT  

# Analyze and Visualize the EPM search results   ############################################  
mo_proj -------------------------  
After searching for motif matches within a reference sequence, these results can be tested and further characterized using the occ_filter, mo_feat-filter, mo_feature_tester and mo_predictability scripts.  
#When the occurence files are to large please use the split_files.sh to split them #########  

./mo_proj/occ_filter_v1.1R  
./mo_proj/mo_feat-filter.v3.4.R  
./mo_proj/mo_feature_tester.v1.0.R  
./mo_proj/mo_predictabilityV1.5.R  
./mo_proj/mo_check_mapping-performance_V1.7.R  
./mo_proj/mo_genotype_variance.v1.4.R  
  
##############################################################################################  
ref_seq -------------------------  
Store genome data like fasta, gff and many more here for INPUT.   
0MODELS (formerly 0MOTIFS)  
Store deepCRE output here.   
##############################################################################################  



---------- Overview --------  
recommended directory structure  
Please find more detailed descriptions  
of the directories and the code in the files itself.  

This is a pipeline of consecutive operations that can be and will be availabe here.  

INPUT DIRECTORY    /ref_seq                 /0MODELS                               
                 - fastas [deepCRE input]   - HDF5/H5.file [deepCRE output]         
                 - gffs           _____________|________________________________                       
                 - meta-data      |                      |                     |  
START DIRECTORY      |        /mo_nom                   /mo_range              /mo_imp  
output               |      - get motif patterns      - get motif meta-data   - get model meta data  
                     |      - motif annotation           |                     |    
                     |      - motif modification         |                     |     
                     |____________|______________________|_____________________|  
                                  |                      |  
                              /mo_clu                    |        MAPPING to reference (external)  
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








###################################################################################################

mo_old ---------------------------
Old and outdated scripts used for the moca_blue suite are stored here.






