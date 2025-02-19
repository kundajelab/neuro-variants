# SLURM command to request resources for an interactive bash session
# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# Activate the conda environment containing R:
# conda activate r
# Start R:
# R

################################################################################

# load
library(data.table)

# input arguments
variantSet="rare"
bias="K562_bias"

# for (variantSet in c("common","rare")) {
for (variantSet in c("asd")) {
  
  # PhyloP:
  if (variantSet=="asd") {
    f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/Trevino_et_al_AllMutations.chrALL.phylop.bed"
  } else if (variantSet=="rare") {
    f="/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chrALL.filter.score.v2.phylop.bed"
  } else if (variantSet=="common") {
    f="/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.v2.phylop.bed"
    #need to do
  }
  phylop.df = fread(f,data.table = F,stringsAsFactors = F)
  colnames(phylop.df) = c("chr",'pos0',"snp_pos","snp_id","phylop")
  phylop.df$snp_id2 = unlist(lapply(strsplit(phylop.df$snp_id,":"),function(x) substring(paste(x[c(1,2,4)],collapse = "_"),4)))
  
  # CADD:
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/annotations/CADD.",variantSet,".non_coding.tsv")
  cadd = fread(f,data.table = F,stringsAsFactors = F)
  cadd = cadd[,!(colnames(cadd) %in% c("POS","REF","ALT","CHROM","SubjectID","LoF"))]
  
  # chromBPnet
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".all_dataset.K562_bias.txt")
  cbp.all = fread(f,data.table = F,stringsAsFactors = F)
  # idx = anyDuplicated(colnames(cbp.all)); if (idx != 0) {cbp.all = cbp.all[,-idx]}
  
  df = merge(phylop.df[,c("snp_id","snp_id2","phylop")],cadd,by.x=c("snp_id2"),by.y=c("id"))
  df = merge(df,cbp.all,by.x=c("snp_id"),by.y=c("variant_id"))

  ##############################################################################
  
  # Read gene constraint scores:
  
  # pLOF (gnomad): (lower LOUEF values is more constrained)
  plof = fread("/oak/stanford/groups/smontgom/amarder/bin/gnomad.v2.1.1.lof_metrics.by_gene.txt",data.table = F,stringsAsFactors = F)
  plof = plof[,c("gene","oe_lof_upper","oe_lof_upper_bin","pLI")]
  
  f="/oak/stanford/groups/smontgom/amarder/bin/shet_zeng_et_al.csv"
  shet = fread(f,data.table = F,stringsAsFactors = F)
  shet = shet[,c('ensg',"post_mean")]
  
  ensg_hgnc = fread("/oak/stanford/groups/smontgom/amarder/data/ensg_hgnc_mapping.v87.txt",data.table = F,stringsAsFactors = F)
  plof.tmp = merge(plof,ensg_hgnc,by.x="gene",by.y="hgnc")
  shet.tmp = merge(shet,ensg_hgnc,by.x="ensg",by.y="ensg")
  plof.tmp$keep = plof.tmp$ensg %in% shet.tmp$ensg
  plof.tmp = plof.tmp[order(-1*plof.tmp$keep,plof.tmp$oe_lof_upper,decreasing = F),]
  shet.tmp$keep = shet.tmp$hgnc %in% plof.tmp$gene
  shet.tmp = shet.tmp[order(shet.tmp$keep,shet.tmp$post_mean,decreasing = T),]
  plof.tmp = plof.tmp[!(duplicated(plof.tmp$gene)),]
  shet.tmp = shet.tmp[!(duplicated(shet.tmp$hgnc)),]
  # gene_constraint_table = merge(plof.tmp,shet.tmp,by.x="gene",by.y="hgnc",all=T)

  # Only merge shet, ignore plof:
  shet.tmp1 = shet.tmp[,c("hgnc","post_mean")]; colnames(shet.tmp1) = c("hgnc","s_het_1")
  shet.tmp2 = shet.tmp[,c("hgnc","post_mean")]; colnames(shet.tmp2) = c("hgnc","s_het_2")
  shet.tmp3 = shet.tmp[,c("hgnc","post_mean")]; colnames(shet.tmp3) = c("hgnc","s_het_3")
  df = merge(df,shet.tmp1,by.x="closest_gene_1",by.y="hgnc",all.x=TRUE)
  df = merge(df,shet.tmp2,by.x="closest_gene_2",by.y="hgnc",all.x=TRUE)
  df = merge(df,shet.tmp3,by.x="closest_gene_3",by.y="hgnc",all.x=TRUE)
  
  # # compared to matching on ENSG ensembl id from CADD, there is a cor of ~0.7 in s_het values
  # # however, 2x more matching using hgnc symbols instead!
  # tmp = merge(df[,c("snp_id","GeneName","s_het_1","s_het_2","s_het_3")],shet.tmp[,c("ensg","post_mean")],by.x="GeneName",by.y="ensg")
  # tmp = tmp[!duplicated(tmp$snp_id),]
  # cor.test(tmp$s_het_1,tmp$post_mean,use='pairwise.complete.obs')
  # sum(is.na(df$s_het_1))
  # dim(tmp)
  
  # 1 vs 2 or is correlated:
  # cor.test(df$s_het_3,df$s_het_2,use='pairwise.complete.obs')

  # currently ignoring plof
  # df$loeuf = df$oe_lof_upper 
  # df = merge(df,plof[,c("gene","oe_lof_upper","oe_lof_upper_bin","pLI")],by.x="closest_gene_1",by.y="gene",all.x = TRUE)
  
  #######################################################################################
  
  # this is script to remove outliers - need to finish!
  f = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/model_performance.outliers.tsv"

  # Read in outliers file
  outliers <- fread(f, data.table = FALSE, stringsAsFactors = FALSE)
  
  # Generate unique patterns
  cell <- paste0(outliers$celltype, ".", outliers$dataset)
  patterns <- unique(cell)
  
  # Find indices of outliers in df
  indices <- unlist(sapply(patterns, function(p) grep(p, colnames(df))))
  indices <- as.vector(indices) # Ensure it is a simple vector
  
  # Remove outlier columns from df
  df <- df[, -indices]
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot1.txt")
  # fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  
  #######################################################################################
  
  # If running in 2 parts, read in first annotation part:
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot1.txt")
  # df = fread(f,data.table = F,stringsAsFactors = F)
  
  # impute missing gene constraint scores:
  df$s_het_1[is.na(df$s_het_1)] = mean(df$s_het_1,na.rm=TRUE)
  df$s_het_2[is.na(df$s_het_2)] = mean(df$s_het_2,na.rm=TRUE)
  df$s_het_3[is.na(df$s_het_3)] = mean(df$s_het_3,na.rm=TRUE)
  
  # mean s_het
  df$s_het_avg = apply(df[,c("s_het_1","s_het_2","s_het_3")],1,mean)

  # chrombpnet summary info, such as max scores or num cell types affected
  df$cbp_min_pval = apply(df[,grep("abs_logfc.mean.pval",colnames(df))],1,min,na.rm=T)
  cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
  cols <- grep("pval", cols, invert = TRUE, value = TRUE)
  df$cbp_max_score = apply(df[,cols],1,max,na.rm=T)
  df$num_peaks = apply(df[,grep("peak_overlap.",colnames(df))],1,sum,na.rm=T)
  df$num_cbp = apply(df[,grep("abs_logfc.mean.pval",colnames(df))],1,function(x) {sum(x<0.01,na.rm = T)})
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot2.txt")
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
}

