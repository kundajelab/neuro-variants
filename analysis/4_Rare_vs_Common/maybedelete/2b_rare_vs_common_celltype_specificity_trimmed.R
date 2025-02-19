  
# repeat cell type specificity analysis, except use trimmed cell types
context="adult_brain"
adult_brain = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.",context,".txt"),data.table = F,stringsAsFactors = F,header = F)[,1]
context="fetal_brain"
fetal_brain = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.",context,".txt"),data.table = F,stringsAsFactors = F,header = F)[,1]
context="adult_heart"
adult_heart = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.",context,".txt"),data.table = F,stringsAsFactors = F,header = F)[,1]
context="fetal_heart"
fetal_heart = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.",context,".txt"),data.table = F,stringsAsFactors = F,header = F)[,1]

keep = adult_brain
df$num_peaks_ab = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_ab = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = fetal_brain
df$num_peaks_fb = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_fb = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = fetal_heart
df$num_peaks_fh = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
df$num_peakscbp_fh = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})

keep = adult_heart
df$num_peaks_ah = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
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
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/FinalAnalysis/rare_vs_common.celltype_specificity_trimmed.txt")
fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

  
  