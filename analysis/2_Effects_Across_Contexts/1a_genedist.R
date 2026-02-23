# SLURM job submission command for running on HPC cluster
# srun --account=smontgom --partition=batch --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# Activate the R conda environment
# conda activate r
# Launch R
# R

# Load required libraries
library(data.table)  # For fast data manipulation
library(stringr)      # For string operations

# Read a sample of the data (first 10,000 lines)
# df = fread(cmd = "head -10000 /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",
#            data.table = F, stringsAsFactors = F,
#            select = c(
#              paste0("peak_overlap.", keep),
#              paste0("abs_logfc.mean.", keep),
#              paste0("abs_logfc.mean.pval.", keep),
#              "snp_id", "chr", "pos", "gene_distance_1", "s_het_1"
#            )
# )

# Read the full dataset
df = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",
           data.table = F, stringsAsFactors = F)

# Extract column names that contain "peak_overlap."
cell_types = colnames(df)[grep("peak_overlap.", colnames(df))]
# Keep only those containing "domcke"
cell_types = cell_types[grep("domcke", cell_types)]
# Remove "peak_overlap." and ".domcke_2020" suffix
cell_types <- gsub("peak_overlap\\.", "", cell_types)
cell_types <- gsub("\\.domcke_2020", "", cell_types)

# Define the columns to keep based on cell type
keep = paste0(cell_types, ".domcke_2020")

# Extract p-values for peaks
tmp = df[, paste0("abs_logfc.mean.pval.", keep)]
# Set p-value to 1 where there is no peak overlap
for (i in 1:length(keep)) {
  tmp[, i][!df[, paste0("peak_overlap.", keep[i])]] = 1
}

# Count the number of significant cell types that are accessible with chrombpnet p < 0.01
df$num_peakscbp = apply(tmp, 1, function(x) { sum(x < 0.01, na.rm = T) })

# Recalculate num_peakscbp with a stricter threshold (p < 0.005)
df$num_peakscbp005 = apply(tmp, 1, function(x) { sum(x < 0.005, na.rm = T) })

# Recalculate num_peakscbp with a stricter threshold (p < 0.1)
df$num_peakscbp1 = apply(tmp, 1, function(x) { sum(x < 0.1, na.rm = T) })

# Recalculate num_peakscbp with a stricter threshold (p < 0.1)
df$num_peakscbp05 = apply(tmp, 1, function(x) { sum(x < 0.05, na.rm = T) })

# Get the maximum absolute log fold change across cell types
df$max_cbp = apply(df[, paste0("abs_logfc.mean.", keep)], 1, max)

# Count the number of overlapping peaks per variant
df$num_peaks = apply(df[, paste0("peak_overlap.", keep)], 1, sum, na.rm = T)

# Subset data to keep only rows where at least one peak is present
df.sub = df[df[, "num_peaks"] > 0,]

# Log-transform gene distance (adding 1 to avoid log(0))
df.sub$gene_distance_1_log10 = log10(df.sub$gene_distance_1 + 1)

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp >= ceiling(max(df.sub$num_peakscbp) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp == 0] = "Null"
df.sub$set[df.sub$num_peakscbp == 1] = "Specific"

# Convert 'set' to a factor with ordered levels
df.sub$set = factor(df.sub$set, levels = c("Null", "Specific", "Multiple", "Shared"))

# Compute mean max_cbp and gene_distance_1 for each 'set' category
aggregate(max_cbp ~ set, df.sub, mean)
aggregate(gene_distance_1 ~ set, df.sub, mean)

# Fit linear models to assess relationships between variables
summary(lm(gene_distance_1_log10 ~ s_het_1 + set, df.sub))
summary(lm(gene_distance_1_log10 ~ s_het_1 + max_cbp, df.sub))
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp, df.sub))

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp1 >= ceiling(max(df.sub$num_peakscbp1) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp1 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp1 == 1] = "Specific"

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp05 >= ceiling(max(df.sub$num_peakscbp05) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp05 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp05 == 1] = "Specific"

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp005 >= ceiling(max(df.sub$num_peakscbp005) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp005 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp005 == 1] = "Specific"
df.sub$set = factor(df.sub$set,levels=c("Null","Specific","Multiple","Shared"))

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))

# Reset the "set" variable to the original:
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp >= ceiling(max(df.sub$num_peakscbp) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp == 0] = "Null"
df.sub$set[df.sub$num_peakscbp == 1] = "Specific"
df.sub$set = factor(df.sub$set, levels = c("Null", "Specific", "Multiple", "Shared"))

# Compute mean absolute log fold change across peaks
df.sub$mu_cbp = apply(df.sub[, paste0("abs_logfc.mean.", keep)], 1, mean)

# Subset data to include only variants with at least one peak + chrombpnet P < 0.05
df.sub2 = subset(df.sub, num_peakscbp > 0)

# Fit linear model using max_cbp as a predictor
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp, df.sub2))$coef["max_cbp", "Pr(>|t|)"]

# Fit linear model using mu_cbp as a predictor
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + mu_cbp, df.sub2))$coef["mu_cbp", "Pr(>|t|)"]

# Further subset to include only variants where s_het_1 is in the top 10%
df.sub2 = subset(df.sub, s_het_1 > quantile(s_het_1, probs = 0.9))

# Fit linear model for variants near highly constrained genes
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp, df.sub2))$coef["max_cbp", "Pr(>|t|)"]

# Define output file path
f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"

# Save selected columns to file
fwrite(df.sub[, c("snp_id", "chr", "pos", "set", "closest_gene_1","gene_distance_1_log10", "s_het_1", "max_cbp", "num_peakscbp","num_peakscbp1","num_peakscbp05","num_peakscbp005","phylop")],
       f.out, quote = F, na = "NA", sep = '\t', row.names = F, col.names = T)
