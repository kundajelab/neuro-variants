# conda activate r

library(data.table)
library(stringr)

f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$gene_distance_1_log10 = log10(df$gene_distance_1+1)

constraint_vs_development_specificity = function(label,adultCol,fetalCol) {
  df[,label] = df[,paste0("abs_logfc.mean.",adultCol)] - df[,paste0("abs_logfc.mean.",fetalCol)]
  idx = df[,paste0("peak_overlap.",adultCol)] & df[,paste0("peak_overlap.",fetalCol)]
  sum(idx)
  mod <- lm(df$phylop[idx]~scale(df[,label][idx])+df$gene_distance_1_log10[idx]+df$s_het_1[idx])
  return(summary(mod)$coef["scale(df[, label][idx])",])
}

# subset(df,grepl("Excitatory",cell_v2))
# domcke et al. fetal_brain.Excitatory_neurons 
# Cluster1 corces_2020 Isocortical Excitatory
label = "excitatory_neurons_dev"
adultCol = "Cluster1.corces_2020" # Isocortical Excitatory
fetalCol = "fetal_brain.Excitatory_neurons.domcke_2020" # Excitatory neuron

# adultCol = "Cluster3.corces_2020" # Hippocampal Excitatory 1
# fetalCol = "c4.trevino_2021" # Interneuron

constraint_vs_development_specificity(label,adultCol,fetalCol)

df[,label] = df[,paste0("abs_logfc.mean.",adultCol)] - df[,paste0("abs_logfc.mean.",fetalCol)]
idx = df[,paste0("peak_overlap.",adultCol)] & df[,paste0("peak_overlap.",fetalCol)]
idx_a = df[,paste0("abs_logfc.mean.pval.",adultCol)] < 0.01
idx_f = df[,paste0("abs_logfc.mean.pval.",fetalCol)] < 0.01
idx_a_null = df[,paste0("abs_logfc.mean.pval.",adultCol)] > 0.1
idx_f_null = df[,paste0("abs_logfc.mean.pval.",fetalCol)] > 0.1


adult_thres = quantile(df[idx,label],0.9) 
fetal_thres = quantile(df[idx,label],0.1) 

idx_a_spec = (df[,label] > adult_thres) & idx_a & idx_f_null & idx
idx_f_spec = (df[,label] < fetal_thres) & idx_f & idx_a_null & idx

df[idx_f_spec,c("snp_id",paste0("abs_logfc.mean.pval.",adultCol),paste0("abs_logfc.mean.pval.",fetalCol),label)][2,]
sum(idx_a_spec); sum(idx_f_spec)
sum(idx_a & idx); sum(idx_f & idx)

fisher.test(idx_a & idx,idx_f & idx)
fisher.test(idx_a[idx],idx_f[idx])

summary(lm(df$phylop[idx] ~ idx_f_spec[idx] + df$gene_distance_1_log10[idx]+df$s_het_1[idx]))
summary(lm(df$phylop[idx] ~ idx_a_spec[idx] + df$gene_distance_1_log10[idx]+df$s_het_1[idx]))

idx2 = (idx_a_spec | idx_f_spec)
summary(lm(df$phylop[idx2] ~ idx_a_spec[idx2] + df$gene_distance_1_log10[idx2]+df$s_het_1[idx2]))
summary(lm(df$phylop[idx2] ~ idx_f_spec[idx2] + df$gene_distance_1_log10[idx2]+df$s_het_1[idx2]))

df$adult_specific = idx_a_spec
df$fetal_specific = idx_f_spec

tmp = df[idx2,c("snp_id","closest_gene_1","phylop","gene_distance_1_log10","s_het_1","adult_specific","fetal_specific",paste0("abs_logfc.mean.",adultCol),paste0("abs_logfc.mean.",fetalCol),label)]
f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.txt"
fwrite(tmp,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

tmp = df[idx,c("snp_id","closest_gene_1","phylop","gene_distance_1_log10","s_het_1","adult_specific","fetal_specific",paste0("abs_logfc.mean.",adultCol),paste0("abs_logfc.mean.",fetalCol),label)]
f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.accessible.txt"
fwrite(tmp,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

summary(lm(df$phylop[idx2] ~ idx_f_spec[idx2] + df$s_het_1[idx2] + df$gene_distance_1_log10[idx2]))
summary(lm(df$s_het_1[idx2] ~ idx_a_spec[idx2] + df$gene_distance_1_log10[idx2]))





