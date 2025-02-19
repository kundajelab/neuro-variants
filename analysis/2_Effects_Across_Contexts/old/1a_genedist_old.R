# srun --account=smontgom --partition=batch --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# conda activate r
# R

library(data.table)
library(stringr)

df = fread(cmd = "head -10000 /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",
           data.table = F,stringsAsFactors = F,
           select = 
             c(
               paste0("peak_overlap.",keep),paste0("abs_logfc.mean.",keep),paste0("abs_logfc.mean.pval.",keep),"snp_id","chr","pos","gene_distance_1","s_het_1"
             )
)

df = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

cell_types = colnames(df)[grep("peak_overlap.",colnames(df))]
cell_types = cell_types[grep("domcke",cell_types)]
cell_types <- gsub("peak_overlap\\.", "", cell_types)
cell_types <- gsub("\\.domcke_2020", "", cell_types)

keep = paste0(cell_types,".domcke_2020")
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})
df$max_cbp = apply(df[,paste0("abs_logfc.mean.",keep)],1,max)
df$num_peaks = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)

df.sub = df[df[,paste0("num_peaks")]>0,]
df.sub$gene_distance_1_log10 = log10(df.sub$gene_distance_1 + 1)

# annotate
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp >= ceiling(max(df.sub$num_peakscbp)*0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp==0] = "Null"
df.sub$set[df.sub$num_peakscbp==1] = "Specific"
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Restrained","Multiple","Shared"))

aggregate(max_cbp~set,df.sub,mean)
aggregate(gene_distance_1~set,df.sub,mean)

# statistical model:
summary(lm(gene_distance_1_log10 ~ s_het_1 + set,df.sub))
summary(lm(gene_distance_1_log10 ~ s_het_1 + max_cbp,df.sub))
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp,df.sub))

df.sub$mu_cbp = apply(df.sub[,paste0("abs_logfc.mean.",keep)],1,mean)
df.sub2 = subset(df.sub,num_peakscbp > 0)
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp,df.sub2))$coef["max_cbp","Pr(>|t|)"]
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + mu_cbp,df.sub2))$coef["mu_cbp","Pr(>|t|)"]
df.sub2 = subset(df.sub,s_het_1 > quantile(s_het_1,probs=0.9))
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp,df.sub2))["max_cbp","Pr(>|t|)"]

keep = paste0(cell_types,".domcke_2020")
tmp = df.sub[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df.sub[,paste0("peak_overlap.",keep[i])]] = 1}
df.sub$num_peakscbp = apply(tmp,1,function(x) {sum(x<0.001,na.rm = T)})
summary(lm(num_peakscbp ~ s_het_1 + gene_distance_1_log10 + max_cbp,df.sub))

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
fwrite(df.sub[,c("snp_id","chr","pos","set","gene_distance_1_log10","s_het_1","max_cbp","num_peakscbp")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

