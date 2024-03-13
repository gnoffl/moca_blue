#CONVERT .RData to VCF File Format

#############################################################
# Load the data from RData file
load("AthAly.genome.snp.RData")
# Define header for VCF file
vcf_header <- "##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

#############################################################
convert_to_vcf <- function(row, fasta_names) {
    chrom <- as.character(row["chr"])  # Convert to character in case it's a factor
    if (startsWith(chrom, "Chr")) {
        chrom <- substr(chrom, 4, nchar(chrom))  # Remove "Chr" prefix
    }
    pos <- row["Atpos"]  # "Atpos" -> POS
    id <- row["gene"]  # "gene" -> ID
    ref <- row["tha"]  # "tha" -> REF
    alt <- row["lyr"]  # "lyr" -> ALT
    info <- "."  # Get corresponding fasta_name_0 value based on chromosome
    qual <- 100  # QUAL value
    filt <- "PASS"  # FILTER value
    
    return(paste(chrom, pos, id, ref, alt, qual, filt, info, sep = "\t"))
}

##############################################################
# Update fasta names in VCF table
#table_data <- data.frame(
#  fasta_name_0 = c("1 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF", 
#                   "2 dna:chromosome chromosome:TAIR10:2:1:19698289:1 REF", 
#                   "3 dna:chromosome chromosome:TAIR10:3:1:23459830:1 REF", 
#                   "4 dna:chromosome chromosome:TAIR10:4:1:18585056:1 REF", 
#                   "5 dna:chromosome chromosome:TAIR10:5:1:26975502:1 REF"),
#  fasta_name = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
#)
#
#fasta_names <- setNames(table_data$fasta_name_0, table_data$fasta_name)
#
vcf_records <- lapply(seq_len(nrow(AthAly.genome.snp)), function(i) convert_to_vcf(AthAly.genome.snp[i, ], fasta_names))

writeLines(c(vcf_header, unlist(vcf_records)), "AthAly.genome.snp3.vcf")

