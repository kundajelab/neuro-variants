# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# conda activate r
# R

# load
library(data.table)
library(stringr)

# input arguments
variantSet="rare"
bias="K562_bias"

df = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)
df$gene_distance_1_log10 = log10(df$gene_distance_1 + 1)

# extract unique cols (dataset, cell)
string = grep("abs_logfc.mean.pval.",colnames(df),value = TRUE)
dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")

cell_types = colnames(df)[grep("peak_overlap.",colnames(df))]
cell_types <- gsub("peak_overlap\\.", "", cell_types)
ind1 <- grepl("domcke_2020", cell_types) & grepl("fetal_brain", cell_types)
ind2 <- grepl("trevino_2021", cell_types)
fetal_brain <- cell_types[ind1 | ind2]

ind1 <- grepl("corces_2020", cell_types)
adult_brain <- cell_types[ind1]

ind1 <- grepl("domcke_2020", cell_types) & grepl("fetal_brain", cell_types)
ind2 <- grepl("ameen_2022", cell_types)
fetal_heart <- cell_types[ind1 | ind2]

ind1 <- grepl("encode_2024", cell_types)
adult_heart <- cell_types[ind1]

keep = c(adult_brain,adult_heart,fetal_brain,fetal_heart)
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

popu="sas"

i = 0; res.lst = list()
for (pheno in c("num_peakscbp","cbp_max_score")) {
  for (popu in c("afr","amr","asj","eas","fin","nfe","sas")) {
    print(popu)
    mod <- lm(
      # scale(df[,paste0("cbp_max_score")]) ~
      scale(df[,paste0(pheno)]) ~
        ifelse(df[,paste0("gnomADg_AF_joint_",popu)] > 0.01,"Common","Rare") +
        df$gene_distance_1_log10 + 
        df$s_het_1
    )
    print(modcoef <- summary(mod)$coef[2,])
    i = i + 1
    res.lst[[i]] = data.frame(pheno,popu,t(as.data.frame(modcoef)))
  }
}
res = as.data.frame(do.call(rbind,res.lst))

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/population_drift_comparison.txt"
fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


