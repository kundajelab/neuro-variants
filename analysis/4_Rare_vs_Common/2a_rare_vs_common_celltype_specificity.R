# need 256 to merge common and rare:
# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 
# conda activate r

library(data.table)
library(stringr)

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)
rare = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

# Make column names match
common <- common[, colnames(rare), drop = FALSE]
rare <- rare[, colnames(common), drop = FALSE]

common$set = "common"
rare$set = "rare"
# colnames(rare)[colnames(rare)=="cbp_min_pval_a"] = "cbp_min_pval_fb"
df = as.data.frame(rbind(common,rare))
df$gene_distance_1_log10 = log10(df$gene_distance_1 + 1)

# extract unique cols (dataset, cell)
string = grep("abs_logfc.mean.pval.",colnames(df),value = TRUE)
dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")

################################################################################

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

keep = adult_brain
df$num_peaks_ab = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_ab = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_ab = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = fetal_brain
df$num_peaks_fb = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_fb = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_fb = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = fetal_heart
df$num_peaks_fh = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_fh = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_fh = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = adult_heart
df$num_peaks_ah = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_ah = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_ah = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

res.lst = list(); i = 0
mod = (lm(scale(num_peakscbp_fb)~s_het_1 + gene_distance_1_log10 + set,subset(df,num_peakscbp_fb > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Fetal brain",analysis="orig",beta=res[1],se=res[2],pval=res[4])
mod = (lm(scale(num_peakscbp_fb)~s_het_1 + gene_distance_1_log10 + set + num_peaks_fb,subset(df,num_peakscbp_fb > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Fetal brain",analysis="peaks",beta=res[1],se=res[2],pval=res[4])

mod = (lm(scale(num_peakscbp_fh)~s_het_1 + gene_distance_1_log10 + set,subset(df,num_peakscbp_fh > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Fetal heart",analysis="orig",beta=res[1],se=res[2],pval=res[4])
mod = (lm(scale(num_peakscbp_fh)~s_het_1 + gene_distance_1_log10 + set + num_peaks_fh,subset(df,num_peakscbp_fh > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Fetal heart",analysis="peaks",beta=res[1],se=res[2],pval=res[4])

mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + set,subset(df,num_peakscbp_ab > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Adult brain",analysis="orig",beta=res[1],se=res[2],pval=res[4])
mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + set + num_peaks_ab,subset(df,num_peakscbp_ab > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Adult brain",analysis="peaks",beta=res[1],se=res[2],pval=res[4])

mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + set,subset(df,num_peakscbp_ah > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Adult heart",analysis="orig",beta=res[1],se=res[2],pval=res[4])
mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + set + num_peaks_ah,subset(df,num_peakscbp_ah > 0)))
res = summary(mod)$coef["setrare",]
i=i + 1; res.lst[[i]] = data.frame(context="Adult heart",analysis="peaks",beta=res[1],se=res[2],pval=res[4])

res = as.data.frame(do.call(rbind,res.lst))
res

# save results:
f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.celltype_specificity.txt")
fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)




