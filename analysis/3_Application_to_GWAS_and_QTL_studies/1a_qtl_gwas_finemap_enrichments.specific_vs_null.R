library(data.table)
library(stringr)

################################################################################

# Read in common variant data:

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

###################

# Extract column names that contain "peak_overlap."
cell_types = colnames(common)[grep("peak_overlap.", colnames(common))]
# Keep only those containing "domcke"
cell_types = cell_types[grep("domcke", cell_types)]
cell_types = cell_types[grep("fetal_brain", cell_types)]
# Remove "peak_overlap." and ".domcke_2020" suffix
cell_types <- gsub("peak_overlap\\.", "", cell_types)
cell_types <- gsub("\\.domcke_2020", "", cell_types)

# Define the columns to keep based on cell type
keep = paste0(cell_types, ".domcke_2020")
# keep = paste0(cell_types, ".corces_2020")

# Extract p-values for peaks
tmp = common[, paste0("abs_logfc.mean.pval.", keep)]
# Set p-value to 1 where there is no peak overlap
for (i in 1:length(keep)) {
  tmp[, i][!common[, paste0("peak_overlap.", keep[i])]] = 1
}

# Count the number of significant cell types that are accessible with chrombpnet p < 0.01
common$num_peakscbp = apply(tmp, 1, function(x) { sum(x < 0.01, na.rm = T) })

# Count the number of overlapping peaks per variant
common$num_peaks = apply(common[, paste0("peak_overlap.", keep)], 1, sum, na.rm = T)

# # Subset data to keep only rows where at least one peak is present
# common.sub = common[common[, "num_peaks"] > 0,]

# Annotate each row based on the number of significant peaks
common$set = "Multiple"
common$set[common$num_peakscbp >= ceiling(max(common$num_peakscbp) * 0.8)] = "Shared"
common$set[common$num_peakscbp == 0] = "Null"
common$set[common$num_peakscbp == 1] = "Specific"






################################################################################

# This section analyzes fine-mapped GWAS meta-analysis results that have been fine-mapped using an out-of-sample 1000G reference panel.

organ="brain"
common.sub = subset(common,set %in% c("Null","Specific"))
common.sub = subset(common,set %in% c("Null","Specific") & num_peaks > 0)

finemap = fread(paste0("/oak/stanford/groups/smontgom/amarder/data/finemap_data/gtex/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.",organ,".txt.gz"),data.table = F,stringsAsFactors = F)
finemap$chr_pos <- sub("(_[A-Z]+_[A-Z]+_b38)", "", finemap$V5)
finemap = finemap[order(finemap$V6,decreasing = T),]
finemap = finemap[!duplicated(finemap$chr_pos),]
common.sub$chr_pos = paste0(common.sub$chr,"_",common.sub$pos)

common.sub$gwas = NA
finemap.sub = subset(finemap,V6 > 0.9)
idx = common.sub$chr_pos %in% finemap.sub$chr_pos
common.sub$gwas[idx] = TRUE
finemap.sub = subset(finemap,V6 < 0.01)
idx = common.sub$chr_pos %in% finemap.sub$chr_pos
common.sub$gwas[idx] = FALSE
table(common.sub$set)
table(common.sub$gwas)

idx = !is.na(common.sub$gwas)
# test predicting fine map variants
mod = glm((common.sub[idx,"set"]=="Specific")~common.sub$gwas[idx] + log10(common.sub$gene_distance_1[idx]+1) + common.sub$s_het_1[idx],family = binomial(link="logit"))
res = summary(mod)$coef[2,]
res


