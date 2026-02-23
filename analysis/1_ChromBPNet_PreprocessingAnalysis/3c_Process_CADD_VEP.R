# SLURM command to request resources for an interactive bash session
# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# Activate the conda environment containing R
# conda activate r
# Start R
# R

################################################################################

# Variant Set Loop: Iterates over specified variant sets (asd, rare, common), processing each set's corresponding input file.
# Variant Filtering: Separates coding from non-coding variants based on predefined coding_terms.
# Data Cleaning: Handles missing values, removes duplicates, and converts columns with comma-separated values into numeric format.
# Minor Allele Frequency (MAF): Converts allele frequencies into MAF for further analysis.
# Output: Saves processed coding and non-coding variants to separate output files for downstream applications.

################################################################################

# Load required libraries for data manipulation and string operations
library(data.table)
library(dplyr)
library(stringr)

# Loop through different variant sets to process their respective files
# for (variantSet in c("asd","rare", "common","chd")) {
# for (variantSet in c("ldsc")) {
for (variantSet in c("rosmap")) {
    print(variantSet)  # Print the current variant set being processed
  
  # Assign input file paths based on the variant set
  if (variantSet == "asd") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/Trevino_et_al_AllMutations.sort.all.tsv")
  } else if (variantSet == "rare") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/andrew/chrALL.filter.score_input.", variantSet, ".CADD.VEP.gnomad.split.tsv")
  } else if (variantSet == "common") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/1kg.common.gt_0.05.sort.all.tsv")
  } else if (variantSet == "chd") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/chd_snv_list.sort.all.tsv")
  } else if (variantSet == "ldsc") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/ldsc.filtered.variants.hg38.sort.all.tsv")
  } else if (variantSet == "rosmap") {
    f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/rosmap.sort.all.tsv")
  }
  
  print("Reading data...")  # Indicate data loading step
  df = fread(f, data.table = F, stringsAsFactors = F)  # Load data into a data frame
  print("Data loaded.")  # Confirm data loading
  
  # Keep a full copy of the data
  df.full = df
  
  # Create a unique identifier for each variant based on chromosome, position, and alternate allele
  df$id = paste0(df$CHROM, "_", df$POS, "_", df$ALT)
  
  # Remove redundant columns related to SIFT and PolyPhen annotations
  df = df[, !(colnames(df) %in% c("SIFTval", "SIFTcat", "PolyPhenCat", "PolyPhenVal"))]
  
  # Define a list of coding-related terms
  coding_terms <- c(
    'coding_sequence_variant', 'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost',
    'frameshift_variant', 'inframe_insertion', 'inframe_deletion', 
    'splice_donor_variant', 'splice_acceptor_variant',
    'start_lost',
    "non_coding_transcript_exon_variant",
    "mature_miRNA_variant",
    "start_retained_variant","stop_retained_variant",
    "incomplete_terminal_codon_variant",
    "protein_altering_variant"
  )
  
  # Convert "." (missing values) in gnomAD frequency columns to 0
  gnomad_cols <- grep("gnomAD", colnames(df))  # Identify gnomAD columns
  df[, gnomad_cols] <- lapply(df[, gnomad_cols], function(x) as.numeric(ifelse(x == ".", 0, x)))
  
  # Separate coding variants based on the consequence annotation
  df.coding <- df %>%
    filter(str_detect(Consequence, paste(coding_terms, collapse = "|")))
  df = subset(df, !(id %in% df.coding$id))  # Exclude coding variants
  
  # Load protein-coding gene annotations
  protein_coding_genes = fread("/oak/stanford/groups/smontgom/amarder/data/protein_coding_genes.txt", data.table = F, stringsAsFactors = F)
  
  # Add a column indicating whether the gene is protein-coding
  df$protein_coding = df$GeneName %in% protein_coding_genes$`Gene stable ID`
  
  # Prioritize protein-coding genes and remove duplicate variants
  df = df[order(df$protein_coding, decreasing = T),]
  df = df[!duplicated(df$id),]
  df.coding = df.coding[!duplicated(df.coding$id),]
  
  # Print summary of coding variants removed and non-coding variants retained
  paste0("# coding variants removed: ", nrow(df.coding))
  paste0("# non-coding variants retained: ", nrow(df))
  
  # Clean numerical columns with multiple values separated by commas
  useCol = c("GC", "CpG", "bStatistic", "priPhCons", "mamPhCons", "verPhCons", "priPhyloP", "mamPhyloP", "verPhyloP", "GerpN", "GerpS", "PHRED")
  df[, useCol] = df[, useCol] %>% mutate(across(all_of(useCol), ~ as.numeric(gsub(",.*", "", .))))
  df.coding[, useCol] = df.coding[, useCol] %>% mutate(across(all_of(useCol), ~ as.numeric(gsub(",.*", "", .))))
  
  # Convert gnomAD frequencies to Minor Allele Frequencies (MAF)
  maf_converter = function(values) {
    one_minus_values <- 1 - values
    min_values <- pmin(values, one_minus_values)
    return(min_values)
  }
  df[, grepl("gnomAD", colnames(df))] = apply(df[, grepl("gnomAD", colnames(df))], 2, maf_converter)
  df.coding[, grepl("gnomAD", colnames(df.coding))] = apply(df.coding[, grepl("gnomAD", colnames(df.coding))], 2, maf_converter)
  
  # Save non-coding and coding variant data to output files
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/annotations/CADD.", variantSet, ".non_coding.tsv")
  fwrite(df, f.out, quote = F, na = "NA", sep = '\t', row.names = F, col.names = T)
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/annotations/CADD.", variantSet, ".coding.tsv")
  fwrite(df.coding, f.out, quote = F, na = "NA", sep = '\t', row.names = F, col.names = T)
}
