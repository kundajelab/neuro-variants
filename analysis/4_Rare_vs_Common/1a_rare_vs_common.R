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


res.lst=list();j=1
for (i in 1:length(cell_lst)) {
  
  cell=cell_lst[i]
  dataset=dataset_lst[i]
  print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
  
  idx = as.logical(df[,paste0("peak_overlap.",cell,".",dataset)])
  
  y = scale(df[idx,paste0("abs_logfc.mean.",cell,".",dataset)])
  dist = df$gene_distance_1_log10[idx]
  constraint = df[idx,"s_het_1"]
  set = df[idx,"set"]
  
  res = summary(lm(y~set + constraint + dist))
  res = as.data.frame(t(res$coef["setrare",]))
  
  res$cell=cell
  res$dataset=dataset
  res.lst[[j]] = res
  j = j + 1
}

res.df = as.data.frame(do.call(rbind,res.lst))
res.df = res.df[,c(5:ncol(res.df),1:4)]
res.df[order(res.df[,"Estimate"],decreasing=T),][1,]
res.df[order(res.df[,"Estimate"],decreasing=T),]

f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

################################################################################

res.lst=list();j=1
for (i in 1:length(cell_lst)) {
  
  cell=cell_lst[i]
  dataset=dataset_lst[i]
  print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
  
  peak = (df[,paste0("peak_overlap.",cell,".",dataset)])
  dist = df$gene_distance_1_log10
  constraint = df[,"s_het_1"]
  set = df[,"set"]
  
  res <- summary(glm(peak ~ set + constraint + dist, family = binomial(link = "logit")))

  res = as.data.frame(t(res$coef["setrare",]))
  
  res$cell=cell
  res$dataset=dataset
  res.lst[[j]] = res
  j = j + 1
}

res.df = as.data.frame(do.call(rbind,res.lst))
res.df = res.df[,c(5:ncol(res.df),1:4)]
res.df[order(res.df[,"Estimate"],decreasing=T),]

f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.peaks.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


################################################################################

res.lst=list();j=1
for (i in 1:length(cell_lst)) {
  
  cell=cell_lst[i]
  dataset=dataset_lst[i]
  print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
  
  idx.top = as.logical(df[,paste0("peak_overlap.","fetal_brain.Excitatory_neurons",".","domcke_2020")])
  idx = as.logical(df[,paste0("peak_overlap.",cell,".",dataset)]) & idx.top
  
  y = scale(df[idx,paste0("abs_logfc.mean.",cell,".",dataset)])
  dist = df$gene_distance_1_log10[idx]
  constraint = df[idx,"s_het_1"]
  set = df[idx,"set"]
  
  res = summary(lm(y~set + constraint + dist))
  res = as.data.frame(t(res$coef["setrare",]))
  
  y = scale(df[idx,paste0("abs_logfc.mean.","fetal_brain.Excitatory_neurons",".","domcke_2020")])
  res2 = summary(lm(y~set + constraint + dist))
  res2 = as.data.frame(t(res2$coef["setrare",]))
  diff = res2 - res
  
  est1 = res2[1]
  est2 = res[1]
  se1 = res2[2]
  se2 = res[2]
  
  # Difference in estimates
  diff_est <- est1 - est2
  
  # Standard error of the difference
  diff_se <- sqrt(se1^2 + se2^2)
  
  # Confidence interval (95%)
  z <- qnorm(0.975) # z-value for 95% CI
  ci_lower <- diff_est - z * diff_se
  ci_upper <- diff_est + z * diff_se
  res = data.frame(
    diff = as.numeric(diff_est),
    se = as.numeric(diff_se),
    l = as.numeric(ci_lower),
    h = as.numeric(ci_upper)
  )

  res$cell=cell
  res$dataset=dataset
  res$numpeak = sum(idx)
  res.lst[[j]] = res
  j = j + 1
}

res.df = as.data.frame(do.call(rbind,res.lst))
res.df = res.df[,c(5:ncol(res.df),1:4)]
res.df[order(res.df[,"diff"],decreasing=T),]

f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.diff.fetal_brain.Excitatory_neurons.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


################################################################################
