# SLURM command to request resources for an interactive bash session
# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# Activate the conda environment containing R:
# conda activate r
# Start R:
# R

# I apologize in advance for this script being poorly written
# However I have ran it to replicate some findings, and it does work! 
# (run 10/10/25; originally created in August '25)

################################################################################

# load
library(data.table)

# input arguments
variantSet="common"
bias="K562_bias"

f="/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC.hg38.variant_list.tsv"
variant_lst = fread(f,data.table = F,stringsAsFactors = F)

dflst = list()
for (variantSet in c("common","ldsc")) { 
  print(variantSet)
  # PhyloP:
  if (variantSet=="asd") {
    f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/Trevino_et_al_AllMutations.chrALL.phylop.bed"
  } else if (variantSet=="rare") {
    f="/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chrALL.filter.score.v2.phylop.bed"
  } else if (variantSet=="common") {
    f="/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.v2.phylop.bed"
    #need to do
  } else if (variantSet == "chd") {
    # f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/chd_snv_list.sort.all.tsv")
    f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.chrALL.phylop.bed"
  } else if (variantSet == "ldsc") {
    # f = paste0("/oak/stanford/groups/smontgom/erobb/data/watershed/chd_snv_list.sort.all.tsv")
    f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/ldsc.filtered.variants.hg38.chrALL.phylop.bed"
  } else if (variantSet == "rosmap") {
    f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap_variants.chrALL.phylop.bed"
  }

  phylop.df = fread(f,data.table = F,stringsAsFactors = F)
  colnames(phylop.df) = c("chr",'pos0',"snp_pos","snp_id","phylop")
  phylop.df$snp_id2 = unlist(lapply(strsplit(phylop.df$snp_id,":"),function(x) substring(paste(x[c(1,2,4)],collapse = "_"),4)))
  
  # CADD:
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/annotations/CADD.",variantSet,".non_coding.tsv")
  cadd = fread(f,data.table = F,stringsAsFactors = F)
  cadd = cadd[,!(colnames(cadd) %in% c("POS","REF","ALT","CHROM","SubjectID","LoF"))]
  cadd1 = cadd
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/annotations/CADD.",variantSet,".coding.tsv")
  cadd = fread(f,data.table = F,stringsAsFactors = F)
  cadd = cadd[,!(colnames(cadd) %in% c("POS","REF","ALT","CHROM","SubjectID","LoF"))]
  cadd2 = cadd
  cadd2$protein_coding = NA
  cadd = as.data.frame(rbind(cadd1,cadd2))
  
  # chromBPnet
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".all_dataset.K562_bias.txt")
  cbp.all = fread(f,data.table = F,stringsAsFactors = F)
  # idx = anyDuplicated(colnames(cbp.all)); if (idx != 0) {cbp.all = cbp.all[,-idx]}
  
  if (variantSet=="ldsc" | variantSet=="rosmap") {
    # create merge key as idtmp
    phylop.df$idtmp <- paste0(phylop.df$chr, "_", phylop.df$snp_pos)
    # fix off-by-one in position for idtmp
    cadd$idtmp <- unlist(lapply(strsplit(cadd$id, "_"), function(x) {
      chr <- x[1]
      pos <- as.integer(x[2]) - 1
      paste0("chr", chr, "_", pos)
    }))
    dflst[[variantSet]] = merge(phylop.df[,c("snp_id","snp_id2","phylop","idtmp")],cadd,by.x=c("idtmp"),by.y=c("idtmp"))
  } else {
    dflst[[variantSet]] = merge(phylop.df[,c("snp_id","snp_id2","phylop")],cadd,by.x=c("snp_id2"),by.y=c("id"))
  }
  dflst[[variantSet]] = merge(dflst[[variantSet]],cbp.all,by.x=c("snp_id"),by.y=c("variant_id"))
}

dflst[["common"]]$idtmp = paste0(dflst[["common"]]$chr, "_", dflst[["common"]]$pos)
variant_lst$id <- paste(variant_lst[,1], variant_lst[,2], sep = "_")

# sum(dflst[["common"]]$idtmp %in% variant_lst$id)
# sum(dflst[["ldsc"]]$idtmp %in% variant_lst$id)
# nrow(variant_lst)

tmp = subset(dflst[["common"]],idtmp %in% variant_lst$id)
extra_col <- setdiff(colnames(dflst[["ldsc"]]), colnames(tmp))
dflst[["ldsc"]][[extra_col]] <- NULL
variant_lst2 = as.data.frame(rbind(
  tmp,
  dflst[["ldsc"]]
))

variant_lst$set = "annot"

# 2) detect duplicates on the right; decide how to handle
dups_right <- variant_lst2$idtmp[duplicated(variant_lst2$idtmp)]
length(dups_right)  # if > 0, you're in one-to-many land

# 3) rowwise align variant_lst2 to variant_lst by id
idx <- match(variant_lst$id, variant_lst2$idtmp)  # length == nrow(variant_lst)

# optional: drop the join key from the right side to avoid duplication
right_keep <- setdiff(names(variant_lst2), "idtmp")

variant_lst.mg2 <- cbind(
  variant_lst[, c("id", "set")],                 # keep master columns
  variant_lst2[idx, right_keep, drop = FALSE]    # NA where no match
)

# sanity checks
stopifnot(nrow(variant_lst.mg2) == nrow(variant_lst))
stopifnot(identical(variant_lst.mg2$id, variant_lst$id))

f.out=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/ldsc_intermediate_file.annotation.txt")
fwrite(variant_lst.mg2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
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
variant_lst.mg2 = merge(variant_lst.mg2,shet.tmp1,by.x="closest_gene_1",by.y="hgnc",all.x=TRUE,sort=FALSE)
variant_lst.mg2 = merge(variant_lst.mg2,shet.tmp2,by.x="closest_gene_2",by.y="hgnc",all.x=TRUE,sort=FALSE)
variant_lst.mg2 = merge(variant_lst.mg2,shet.tmp3,by.x="closest_gene_3",by.y="hgnc",all.x=TRUE,sort=FALSE)

# # compared to matching on ENSG ensembl id from CADD, there is a cor of ~0.7 in s_het values
# # however, 2x more matching using hgnc symbols instead!
# tmp = merge(variant_lst.mg2[,c("snp_id","GeneName","s_het_1","s_het_2","s_het_3")],shet.tmp[,c("ensg","post_mean")],by.x="GeneName",by.y="ensg")
# tmp = tmp[!duplicated(tmp$snp_id),]
# cor.test(tmp$s_het_1,tmp$post_mean,use='pairwise.complete.obs')
# sum(is.na(variant_lst.mg2$s_het_1))
# dim(tmp)

# 1 vs 2 or is correlated:
# cor.test(variant_lst.mg2$s_het_3,variant_lst.mg2$s_het_2,use='pairwise.complete.obs')

# currently ignoring plof
# variant_lst.mg2$loeuf = variant_lst.mg2$oe_lof_upper 
# variant_lst.mg2 = merge(variant_lst.mg2,plof[,c("gene","oe_lof_upper","oe_lof_upper_bin","pLI")],by.x="closest_gene_1",by.y="gene",all.x = TRUE)

#######################################################################################

# this is script to remove outliers - need to finish!
f = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/model_performance.outliers.tsv"

# Read in outliers file
outliers <- fread(f, data.table = FALSE, stringsAsFactors = FALSE)

# Generate unique patterns
cell <- paste0(outliers$celltype, ".", outliers$dataset)
patterns <- unique(cell)

# Find indices of outliers in variant_lst.mg2
indices <- unlist(sapply(patterns, function(p) grep(p, colnames(variant_lst.mg2))))
indices <- as.vector(indices) # Ensure it is a simple vector

# Remove outlier columns from variant_lst.mg2
variant_lst.mg2 <- variant_lst.mg2[, -indices]

f.out=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot1.txt")
# fwrite(variant_lst.mg2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


#######################################################################################

# If running in 2 parts, read in first annotation part:
f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot1.txt")
# variant_lst.mg2 = fread(f,data.table = F,stringsAsFactors = F)

# impute missing gene constraint scores:
variant_lst.mg2$s_het_1[is.na(variant_lst.mg2$s_het_1)] = mean(variant_lst.mg2$s_het_1,na.rm=TRUE)
variant_lst.mg2$s_het_2[is.na(variant_lst.mg2$s_het_2)] = mean(variant_lst.mg2$s_het_2,na.rm=TRUE)
variant_lst.mg2$s_het_3[is.na(variant_lst.mg2$s_het_3)] = mean(variant_lst.mg2$s_het_3,na.rm=TRUE)

# mean s_het
variant_lst.mg2$s_het_avg = apply(variant_lst.mg2[,c("s_het_1","s_het_2","s_het_3")],1,mean)

# chrombpnet summary info, such as max scores or num cell types affected
variant_lst.mg2$cbp_min_pval = apply(variant_lst.mg2[,grep("abs_logfc.mean.pval",colnames(variant_lst.mg2))],1,min,na.rm=T)
cols <- grep("abs_logfc.mean", colnames(variant_lst.mg2), value = TRUE)
cols <- grep("pval", cols, invert = TRUE, value = TRUE)
variant_lst.mg2$cbp_max_score = apply(variant_lst.mg2[,cols],1,max,na.rm=T)
variant_lst.mg2$num_peaks = apply(variant_lst.mg2[,grep("peak_overlap.",colnames(variant_lst.mg2))],1,sum,na.rm=T)
variant_lst.mg2$num_cbp = apply(variant_lst.mg2[,grep("abs_logfc.mean.pval",colnames(variant_lst.mg2))],1,function(x) {sum(x<0.01,na.rm = T)})

f.out=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot2.txt")
fwrite(variant_lst.mg2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


####################################

# build ldsc annotations!
# for (annot_name in c("plus_cbp","neg_cbp")) {
#   
#   if (annot_name=="plus_cbp") {
#     tmp = data.frame(
#       id=variant_lst.mg2$id,
#       ANNOT=variant_lst.mg2$abs_logfc.mean.Cluster24.corces_2020 < 0.01 & variant_lst.mg2$peak_overlap.Cluster24.corces_2020
#     )
#   } else if (annot_name=="neg_cbp") {
#     tmp = data.frame(
#       id=variant_lst.mg2$id,
#       ANNOT=variant_lst.mg2$abs_logfc.mean.Cluster24.corces_2020 > 0.01 & variant_lst.mg2$peak_overlap.Cluster24.corces_2020
#     )
#   }
#   
#   for (chrnum in 1:22) {
#     print(chrnum)
# 
#     df1 = fread(paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_bed/chr",chrnum,".bed"),data.table = F,stringsAsFactors = F)
#     df2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_lifted/chr",chrnum,".lifted.bed"),data.table = F,stringsAsFactors = F)
#     
#     df2$chrpos = paste0(df2$V1,"_",df2$V3)
#     
#     # 3) rowwise align variant_lst2 to variant_lst by id
#     idx <- match(df2$chrpos, tmp$id)  # length == nrow(variant_lst)
#     
#     # optional: drop the join key from the right side to avoid duplication
#     right_keep <- setdiff(names(tmp), "id")
#     
#     df2b <- cbind(
#       df2,                 # keep master columns
#       tmp[idx, right_keep, drop = FALSE]    # NA where no match
#     )
#     df2b$annot[is.na(df2b$annot)] = FALSE
#     colnames(df2b)[colnames(df2b)=="annot"] = "ANNOT"
#     
#     # sanity checks
#     stopifnot(nrow(df2) == nrow(tmp[idx,]))
#     stopifnot(identical(tmp$id[idx], df2$chrpos))
#     
#     df2b$annot = as.numeric(df2b$annot)
#     
#     idx <- match(df1$V4, df2b$V4)
#     right_keep <- setdiff(names(df2b), "id")
#     df1b <- cbind(
#       df1,                 # keep master columns
#       df2b[idx, right_keep, drop = FALSE]    # NA where no match
#     )
#     
#     f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended/",annot_name,".chr",chrnum,".annot.gz")
#     fwrite(data.frame(ANNOT=df1b$ANNOT),
#            f.out,
#            quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
#   }
# }



#############################################################################################

# SLURM command to request resources for an interactive bash session
# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# Activate the conda environment containing R:
# conda activate r
# Start R:
# R

library(data.table)
f=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/ldsc.all_dataset.K562_bias.annot2.txt")
variant_lst.mg2=fread(f,data.table = F,stringsAsFactors = F)
for (annot_name in c("peak")) {
# for (annot_name in c("plus_cbp","neg_cbp","peak")) {
  if (annot_name == "plus_cbp") {
    tmp <- data.frame(
      # id    = variant_lst.mg2$id,
      # ANNOT = (variant_lst.mg2$abs_logfc.mean.Cluster24.corces_2020 < 0.01) &
      #   (variant_lst.mg2$peak_overlap.Cluster24.corces_2020)
      id    = variant_lst.mg2$id,
      ANNOT = (variant_lst.mg2$abs_logfc.mean.c5.trevino_2021 < 0.01) &
        (variant_lst.mg2$peak_overlap.c5.trevino_2021)
    )
  } else if (annot_name == "neg_cbp") {
    tmp <- data.frame(
      # id    = variant_lst.mg2$id,
      # ANNOT = (variant_lst.mg2$abs_logfc.mean.Cluster24.corces_2020 > 0.01) &
      #   (variant_lst.mg2$peak_overlap.Cluster24.corces_2020)
      id    = variant_lst.mg2$id,
      ANNOT = (variant_lst.mg2$abs_logfc.mean.c5.trevino_2021 > 0.01) &
        (variant_lst.mg2$peak_overlap.c5.trevino_2021)
    )
  } else {
    # tmp <- data.frame(
    #   id    = variant_lst.mg2$id,
    #   ANNOT = (variant_lst.mg2$peak_overlap.Cluster24.corces_2020)
    # )
    tmp <- data.frame(
      id    = variant_lst.mg2$id,
      ANNOT = (variant_lst.mg2$peak_overlap.c5.trevino_2021)
    )
  }
  
  for (chrnum in 1:22) {
    message("chr ", chrnum)
    
    df1 <- fread(paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_bed/chr", chrnum, ".bed"),
                 data.table = FALSE, stringsAsFactors = FALSE)
    df2 <- fread(paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_lifted/chr", chrnum, ".lifted.bed"),
                 data.table = FALSE, stringsAsFactors = FALSE)
    
    df2$chrpos <- paste0(df2$V1, "_", df2$V3)
    
    ## align tmp (by id) onto df2
    idx <- match(df2$chrpos, tmp$id)  # length == nrow(df2)
    right_keep <- setdiff(names(tmp), "id")
    df2b <- cbind(df2, tmp[idx, right_keep, drop = FALSE])
    
    ## fill NAs to 0 and coerce to numeric 0/1
    df2b$ANNOT[is.na(df2b$ANNOT)] <- 0L
    df2b$ANNOT <- as.integer(df2b$ANNOT)
    
    ## sanity checks
    stopifnot(nrow(df2) == length(idx))
    stopifnot(identical(tmp$id[idx], df2$chrpos))  # NAs ok if in same places
    
    ## align df2b onto df1 by SNP ID (V4)
    m <- match(df1$V4, df2b$V4)
    df1b <- cbind(df1, df2b[m, c("ANNOT"), drop = FALSE])
    df1b$ANNOT[is.na(df1b$ANNOT)] <- 0L
    
    ## write single-column annot (rows aligned to df1/.bim order)
    # f.out <- paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended/",
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended2/",
                    annot_name, ".chr", chrnum, ".annot.gz")
    fwrite(data.frame(ANNOT = df1b$ANNOT), f.out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/","ldsc",".FLARE.2.txt")
flare = fread(f,data.table = F,stringsAsFactors = F)



# for (annot_name in c("plus_cbp","neg_cbp")) {
#   for (chrnum in 1:22) {
#     message("chr ", chrnum)
#     
# echo "running ldsc..."
# path_to_ldsc=/oak/stanford/groups/smontgom/amarder/bin/ldsc
# ldsc=$path_to_ldsc/ldsc.py
# hapmap=$path_to_ldsc/hapmap/hapmap3_snps/hm.${chrnum}.snp
# geno="/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrnum}"
# annotation_output <- paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended/",
#                 annot_name, ".chr", chrnum, ".annot.gz")
# 
# $ldsc --l2 --bfile $geno --ld-wind-cm 1 --annot $annotation_output.annot.gz --thin-annot --out $annotation_output --print-snps $hapmap
# 


# ## --- paths you likely keep constant ---
# path_to_ldsc <- "/oak/stanford/groups/smontgom/amarder/bin/ldsc"
# ldsc_py      <- file.path(path_to_ldsc, "ldsc.py")
# python_bin   <- "python"  # or your env: "/path/to/conda/envs/ldsc/bin/python"
# 
# hapmap_dir   <- file.path(path_to_ldsc, "hapmap", "hapmap3_snps")
# geno_prefix  <- "/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC"
# annot_dir    <- "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended"
# 
# ## helper to run a command and stop on error
# run_cmd <- function(cmd, args) {
#   message("Running: ", cmd, " ", paste(args, collapse = " "))
#   rc <- system2(cmd, args)
#   if (rc != 0) stop("Command failed with exit code ", rc)
# }
# 
# ## main loop
# for (annot_name in c("plus_cbp", "neg_cbp")) {
#   for (chrnum in 1:22) {
#     message("\n==== ", annot_name, " | chr ", chrnum, " ====")
#     
#     # Inputs
#     hapmap <- file.path(hapmap_dir, paste0("hm.", chrnum, ".snp"))
#     bfile  <- paste0(geno_prefix, ".", chrnum)
#     annot_path <- file.path(annot_dir, sprintf("%s.chr%d.annot.gz", annot_name, chrnum))
#     
#     # Quick checks
#     if (!file.exists(annot_path)) stop("Missing annot: ", annot_path)
#     if (!file.exists(paste0(bfile, ".bed"))) stop("Missing plink bed: ", bfile, ".bed")
#     if (!file.exists(paste0(bfile, ".bim"))) stop("Missing plink bim: ", bfile, ".bim")
#     if (!file.exists(paste0(bfile, ".fam"))) stop("Missing plink fam: ", bfile, ".fam")
#     if (!file.exists(hapmap)) stop("Missing HapMap3 list: ", hapmap)
#     
#     # LDSC output prefix (NO .annot.gz here)
#     out_prefix <- file.path(annot_dir, sprintf("%s.chr%d", annot_name, chrnum))
#     
#     # Build LD scores
#     args <- c("conda activate ldsc; ",
#               ldsc_py, "--l2",
#               "--bfile", bfile,
#               "--ld-wind-cm", "1",
#               "--annot", annot_path,
#               "--thin-annot",
#               "--out", out_prefix,
#               "--print-snps", hapmap
#     )
#     run_cmd(python_bin, args)
#     
#     # Optional: quick confirm of expected outputs
#     expected <- paste0(out_prefix, c(".l2.ldscore.gz", ".l2.M.gz", ".l2.log"))
#     missing  <- expected[!file.exists(expected)]
#     if (length(missing)) warning("Missing expected outputs: ", paste(missing, collapse=", "))
#     message("Done: ", out_prefix)
#   }
# }
# 
