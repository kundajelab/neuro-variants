# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# conda activate r
# R

# load
library(data.table)
library(stringr)

f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
model_meta = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/model_performance.tsv",data.table = F,stringsAsFactors = F)

adult_brain_performance = fread("/oak/stanford/groups/akundaje/projects/neuro-variants/tables/variant_intersect/rare/corces_2020.rare_variants.in_peaks.phylop.obs.pred.tsv",data.table = F,stringsAsFactors = F)
fetal_brain_performance = fread("/oak/stanford/groups/akundaje/projects/neuro-variants/tables/variant_intersect/rare/trevino_2021.rare_variants.in_peaks.phylop.obs.pred.tsv",data.table = F,stringsAsFactors = F)

# extract unique cols (dataset, cell)
string = grep("abs_logfc.mean.pval.",colnames(df),value = TRUE)
dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")

evo_metric="phylop"
cell="c4"
dataset="trevino_2021"
# cell="Cluster11"
# dataset="corces_2020"

idx = as.logical(df[,paste0("peak_overlap.",cell,".",dataset)])
df.sub = df[idx,]
performance.sub = subset(fetal_brain_performance,sample==paste0(dataset,".",cell))
# performance.sub = subset(adult_brain_performance,sample==paste0(dataset,".",cell))
setDT(performance.sub)
performance.sub[, phylop := as.numeric(phylop)]
performance.sub$phylop[is.na(performance.sub$phylop)] <- mean(performance.sub$phylop, na.rm = TRUE)
performance.sub <- performance.sub[, 
                                                           .(predicted = mean(predicted, na.rm = TRUE), 
                                                             observed  = mean(observed,  na.rm = TRUE), 
                                                             phylop    = mean(phylop,    na.rm = TRUE)), 
                                                           by = variant_id]
colnames(performance.sub)[colnames(performance.sub)=="phylop"] = "phylop_peak"

df.sub.mg = merge(df.sub,performance.sub,by.x="snp_id",by.y="variant_id")

df.sub.mg$error_abs <- abs(
  exp(df.sub.mg$observed) - exp(df.sub.mg$predicted)
)

cor.test(df.sub.mg$error_abs,df.sub.mg$phylop_peak)

error_abs = scale(df.sub.mg$error_abs)
obs_ct = scale(df.sub.mg$observed)
# pred_ct = scale(df.sub.mg$predicted)
phylop_peak = df.sub.mg$phylop_peak
x = scale(df.sub.mg[,paste0("abs_logfc.mean.",cell,".",dataset)])
phylop = (df.sub.mg[,evo_metric])
dist = log10(df.sub.mg[,"gene_distance_1"]+1)
pval = df.sub.mg[,paste0("abs_logfc.mean.pval.",cell,".",dataset)]
constraint = df.sub.mg[,"s_het_1"]

res=summary(lm(phylop~x+dist+constraint + error_abs + obs_ct + phylop_peak));res
res=summary(lm(phylop~x+dist+constraint));res
res=summary(lm(phylop~x+dist+constraint + error_abs));res

res=summary(lm(phylop~x+dist+constraint + error_abs + obs_ct + phylop_peak));res
res=summary(lm(phylop~x+dist+constraint + error_abs));res
res = as.data.frame(t(res$coef["x",]));res


res.lst=list();j=1; fisher.df = list()
for (i in 1:length(cell_lst)) {
  cell=cell_lst[i]
  dataset=dataset_lst[i]
  print(paste0(i,"/",length(cell_lst),": ",cell))
  
  idx = as.logical(df[,paste0("peak_overlap.",cell,".",dataset)])
  print(paste0("evo: ",evo_metric))
  
  x = scale(df[idx,paste0("abs_logfc.mean.",cell,".",dataset)])
  phylop = (df[idx,evo_metric])
  dist = log10(df[idx,"gene_distance_1"]+1)
  pval = df[idx,paste0("abs_logfc.mean.pval.",cell,".",dataset)]
  constraint = df[idx,"s_het_1"]
  
  for (modelNum in c(5)) {
    if (modelNum==1) {
      res=summary(lm(phylop~x))
      res = as.data.frame(t(res$coef["x",]))
    } else if (modelNum==2) {
      res=summary(lm(phylop~x+dist))
      res = as.data.frame(t(res$coef["x",]))
    } else if (modelNum==3) {
      res=summary(lm(x~dist+phylop))
      res = as.data.frame(t(res$coef["dist",]))
    } else if (modelNum==4) {
      res=summary(lm(phylop~dist))
      res = as.data.frame(t(res$coef["dist",]))
    } else if (modelNum==5) {
      res=summary(lm(phylop~x+dist+constraint))
      res = as.data.frame(t(res$coef["x",]))
    } 
    
    res$model = modelNum
    res$evo_metric = evo_metric
    res$cell=cell
    res$dataset=dataset
    res.lst[[j]] = res
    j = j + 1
    
    # k = 0; fisher.lst = list();
    # for (phylop_threshold in seq(0,4,by=0.25)) {
    #   for (cbp_log10threshold in seq(-1,-5,by=-0.25)) {
    #     k = k + 1
    #     # print(k)
    #     x1 = phylop > phylop_threshold
    #     x2 = pval < 10^(cbp_log10threshold)
    #     if ((length(table(x1))==2) & length(table(x2))==2) { 
    #       fisher_res = fisher.test(x1,x2)
    #       fisher.lst[[k]] = data.frame(cell,
    #                                    phylop_threshold,
    #                                    cbp_log10threshold,
    #                                    or=as.numeric(fisher_res$estimate),
    #                                    l=fisher_res$conf.int[1],
    #                                    h=fisher_res$conf.int[2],
    #                                    pval=fisher_res$p.value)
    #     }
    #   }
    # }
    # fisher.df[[cell]] = as.data.frame(do.call(rbind,fisher.lst))
    
  }
}
res.df = as.data.frame(do.call(rbind,res.lst))
res.df = res.df[,c(5:ncol(res.df),1:4)]
# fisher.df.all = as.data.frame(do.call(rbind,fisher.df))

res.df[order(res.df$Estimate)[1:3],]

f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/phylop_vs_cbp.lm.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/phylop_vs_cbp.fisher.txt")
# fwrite(fisher.df.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# 
# res.df.sub = subset(res.df,dataset%in%c("corces_2020","trevino_2021") & model==1)
# t.test(Estimate~dataset,res.df.sub)
# 
# res.df.sub = subset(res.df,dataset%in%c("corces_2020","trevino_2021") & model==2)
# t.test(Estimate~dataset,res.df.sub)
# 
# res.df.sub = subset(res.df,dataset%in%c("corces_2020","trevino_2021") & model==5)
# t.test(Estimate~dataset,res.df.sub)
# 
# # res.df.sub = subset(res.df,dataset%in%c("corces_2020","trevino_2021") & model==4)
# # t.test(Estimate~dataset,res.df.sub)
# 
# f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/phylop_vs_cbp.lm.txt")
# fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# 
# f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/phylop_vs_cbp.fisher.txt")
# fwrite(fisher.df.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
