for chr in {1..22}; do
echo $chr
in_file=/oak/stanford/groups/smontgom/erobb/data/xqtl/watershed/ROSMAP.af.snp.p0.0027.eOutliers.sOutliers.pOutliers.chr${chr}.all.tsv
out_file=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/rosmap.af.snp.p0.0027.eOutliers.sOutliers.pOutliers.chr${chr}.all.trimmed.tsv

cut -f1-6,32-34 "$in_file" > "$out_file"
done

##################

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r


library(data.table)

# ---- Paths ----
pred_path <- "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/rosmap.FLARE.outliers.txt"
df_template <- "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/rosmap.af.snp.p0.0027.eOutliers.sOutliers.pOutliers.chr%s.all.trimmed.tsv"
out_dir <- "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load 'pred' once and prepare IDs ----
pred <- fread(pred_path, data.table = TRUE)
# Make a clean 'id' like "22:17080378" from columns chr (e.g., "chr22") and pos
pred[, id := paste0(sub("^chr", "", as.character(chr)), ":", pos)]

# (Optional but handy) set keys for faster merges
setkey(pred, id)

merged_list <- vector("list", 22)

for (chrnum in 1:22) {
  message(sprintf("Processing chr%s ...", chrnum))
  
  # ---- Load per-chromosome df ----
  df_path <- sprintf(df_template, chrnum)
  df <- fread(df_path, data.table = TRUE)
  
  # Robust CHROM handling (works if CHROM is "22" or "chr22")
  df[, CHROM_clean := sub("^chr", "", as.character(CHROM))]
  df[, id := paste0(CHROM_clean, ":", POS)]
  
  # convert signed P to z scores
  signed_p_to_z <- function(p_signed) {
    sign(p_signed) * qnorm(pmin(1, abs(p_signed)/2), lower.tail = FALSE)
  }
  df$z_eOutliers = signed_p_to_z(df$eOutliers)
  df$z_sOutliers = signed_p_to_z(df$sOutliers)
  df$z_pOutliers = signed_p_to_z(df$pOutliers)
  
  # which columns do you want to average? 
  cols_to_avg <- c("z_eOutliers", "z_sOutliers", "z_pOutliers")
  df_avg <- df[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
               by = .(GeneName, id), 
               .SDcols = cols_to_avg]

  # ---- Subset pred to this chromosome and merge ----
  pred_sub = subset(pred,chr==paste0("chr", chrnum))
  # Merge on id; allow.cartesian=TRUE in case of many-to-one joins (usually not needed, but safe)
  df_mg <- merge(df_avg, pred_sub, by = "id", allow.cartesian = TRUE)
  
  # ---- Write per-chromosome merged file ----
  out_chr_path <- file.path(out_dir, sprintf("rosmap.FLARE.outliers.merged.chr%s.tsv", chrnum))
  fwrite(df_mg, out_chr_path, sep = "\t",quote = F,na = "NA",row.names = F,col.names = T)
  
  # save
  merged_list[[chrnum]] <- df_mg
}

# ---- Combine all chromosomes and write single file ----
merged_all <- rbindlist(merged_list, use.names = TRUE, fill = TRUE)
out_all_path <- file.path(out_dir, "rosmap.FLARE.outliers.merged.all_chr.tsv")
fwrite(merged_all, out_all_path, sep = "\t")

message("Done.")
message(sprintf("All-chromosomes merged file: %s", out_all_path))

#################

library(data.table)
df.mg = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged/rosmap.FLARE.outliers.merged.all_chr.tsv",data.table = F,stringsAsFactors = F)
mapping = fread("/oak/stanford/groups/smontgom/amarder/data/ensg_hgnc_mapping_115.txt",data.table = F,stringsAsFactors = F)
dist_end = abs(mapping$`Transcription start site (TSS)` - mapping$`Gene end (bp)`)
dist_start = abs(mapping$`Transcription start site (TSS)` - mapping$`Gene start (bp)`)
dist_min = apply(cbind(dist_end,dist_start),1,min)
mapping = mapping[order(dist_min,decreasing = FALSE),]
mapping = mapping[!duplicated(mapping$`Gene stable ID`),]

df.mg = merge(df.mg,mapping,by.x="GeneName",by.y="Gene stable ID")
df.mg = subset(df.mg,`Gene name`!="")
df.mg = subset(df.mg,id %in% subset(df.mg,duplicated(id))$id)
dist_end = abs(df.mg$pos - df.mg$`Gene end (bp)`)
dist_start = abs(df.mg$pos - df.mg$`Gene start (bp)`)
dist_min = apply(cbind(dist_end,dist_start),1,min)
df.mg$dist_min = dist_min
dist_tss = abs(df.mg$pos - df.mg$`Transcription start site (TSS)`)
df.mg$dist_tss = dist_tss

z_thres=-3
df.mg.sub = subset(df.mg,FLARE_ab > quantile(FLARE_ab,0.99) & z_eOutliers < z_thres)
df.mg.sub = subset(df.mg,id %in% df.mg.sub$id)
df.mg.sub = df.mg.sub[order(df.mg.sub$id),]

# assume your data.frame is df.mg.sub
d <- df.mg.sub

# impacted if z_eOutliers < -1 (and not NA)
d$impacted <- is.finite(d$z_eOutliers) & d$z_eOutliers < z_thres

# distance from variant to the CURRENT gene on the row
# (your example shows this in 'dist_tss'; change here if you prefer another distance)
# d$dist_to_gene <- d$dist_tss
d$dist_to_gene <- d$dist_min

# helper: summarize per variant id
per_id <- function(g) {
  # handle all-NA distances defensively
  if (all(is.na(g$dist_to_gene))) {
    return(data.frame(
      id = g$id[1],
      closest_gene = NA_character_,
      min_dist = NA_integer_,
      impacted_closest = NA,
      impacted_farther = NA,
      impacted_both = NA,
      impacted_genes = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  min_dist <- min(g$dist_to_gene, na.rm = TRUE)
  
  # genes at the minimal distance (allow ties) 
  is_min <- (!is.na(g$dist_to_gene) & g$dist_to_gene == min_dist) | (g$dist_to_gene < 2000)
  closest_genes <- unique(g[is_min, "Gene name"])
  
  # impacted genes overall
  impacted_genes <- unique(g[g$impacted, "Gene name"])
  
  impacted_closest <- any(g$impacted & is_min, na.rm = TRUE)
  impacted_farther <- any(g$impacted & !is_min, na.rm = TRUE)
  
  data.frame(
    id = g$id[1],
    closest_gene = paste(closest_genes, collapse = ";"),
    min_dist = min_dist,
    impacted_closest = impacted_closest,
    impacted_farther = impacted_farther,
    impacted_both = impacted_closest & impacted_farther,
    impacted_genes = if (length(impacted_genes)) paste(impacted_genes, collapse = ";") else "",
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, lapply(split(d, d$id), per_id))

# 'res' now has one row per variant with flags for closest/farther/both.
# Example peek:
res
table(res$impacted_closest)

tmp = subset(df.mg,id %in% res$id)
tmp[order(tmp$FLARE_ab),]


# subset(df.mg,!(`GeneName` %in% mapping$`Gene stable ID`))[1:5,]

library(data.table)
df.mg = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged/rosmap.FLARE.outliers.merged.all_chr.tsv",data.table = F,stringsAsFactors = F)
setorder(df.mg, +z_eOutliers,na.last = TRUE)
df.mg <- df.mg[!duplicated(df.mg$id),]

rng = seq(0.9,0.999,by=0.001)
zthres.rng = c(-2,-2.5,-3)
iter = 0; res.lst = list()
for (j in 1:length(zthres.rng)) {
  for (i in 1:length(rng)) {
    thres = rng[i]
    zthres = zthres.rng[j]
    cat("FLARE thres: top",thres,"% | outlier z < ",zthres,"\n")
    res = fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,thres),df.mg$z_eOutliers < zthres)
    iter = iter + 1; res.lst[[iter]] = data.frame(iter,thres,zthres,or=res$estimate,l= res$conf.int[1],h=res$conf.int[2],p = res$p.value)
  }
}
res = rbindlist(res.lst)
res
fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/analysis/underexpression_outliers.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

library(data.table)
df.mg = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged/rosmap.FLARE.outliers.merged.all_chr.tsv",data.table = F,stringsAsFactors = F)
setorder(df.mg, -z_eOutliers,na.last = TRUE)
df.mg <- df.mg[!duplicated(df.mg$id),]
rng = seq(0.9,0.999,by=0.001)
zthres.rng = c(3,2.5,2)
iter = 0; res.lst = list()
for (j in 1:length(zthres.rng)) {
  for (i in 1:length(rng)) {
    thres = rng[i]
    zthres = zthres.rng[j]
    cat("FLARE thres: top",thres,"% | outlier z > ",zthres,"\n")
    res = fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,thres),df.mg$z_eOutliers > zthres)
    iter = iter + 1; res.lst[[iter]] = data.frame(iter,thres,zthres,or=res$estimate,l= res$conf.int[1],h=res$conf.int[2],p = res$p.value)
  }
}
res = rbindlist(res.lst)
res
fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/analysis/overexpression_outliers.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

subset(res,zthres==-3)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers > 1)

################################################

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

library(data.table)
f=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/","rosmap",".","all_dataset",".","K562_bias",".annot2.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df.mg = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged/rosmap.FLARE.outliers.merged.all_chr.tsv",data.table = F,stringsAsFactors = F)
mapping = fread("/oak/stanford/groups/smontgom/amarder/data/ensg_hgnc_mapping_115.txt",data.table = F,stringsAsFactors = F)
df.mg.expanded = fread(paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/","rosmap",".FLARE.2.txt"),data.table = F,stringsAsFactors = F)

setorder(df.mg, +z_eOutliers,na.last = TRUE)
df.mg <- df.mg[!duplicated(df.mg$id),]

thres = 0.99
zthres = -3
res = fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,thres),df.mg$z_eOutliers < zthres); res
# table(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,thres),df.mg$z_eOutliers < zthres)

library(dplyr)

# Create an "id" column from snp_id in df.mg.expanded
df.mg.expanded <- df.mg.expanded %>%
  mutate(id = sub("^chr", "", sub(":[ACGT]+:[ACGT]+$", "", snp_id)))

df.mg.expanded = df.mg.expanded[!duplicated(df.mg.expanded$id),]

df <- df %>%
  mutate(id = sub("^chr", "", sub(":[ACGT]+:[ACGT]+$", "", snp_id)))
df = df[!duplicated(df$id),]

# Merge on id
merged <- df.mg %>%
  left_join(df.mg.expanded[,c("id","FLARE_fb_peaks","FLARE_fb")], by = "id")

merged <- merged %>%
  left_join(df[,c("id","peak_overlap.Cluster1.corces_2020","peak_overlap.fetal_brain.Excitatory_neurons.domcke_2020","abs_logfc.mean.Cluster1.corces_2020","abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020")], by = "id")

thres=0.999;zthres = -3
fisher.test(merged$FLARE_ab > quantile(merged$FLARE_ab,thres),merged$z_eOutliers < zthres)

merged$cbp = merged$abs_logfc.mean.Cluster1.corces_2020
merged$cbp_peak = ifelse(merged$peak_overlap.Cluster1.corces_2020, merged$cbp,0)
merged$cbp_nopeak = ifelse(merged$peak_overlap.Cluster1.corces_2020, 0, merged$cbp)
rng = seq(0.9,0.999,by=0.001)
zthres.rng = c(-2,-2.5,-3)
iter = 0; res.lst = list()
for (metric in c("FLARE_ab","cbp","cbp_peak","cbp_nopeak")) {
  for (j in 1:length(zthres.rng)) {
    for (i in 1:length(rng)) {
      thres = rng[i]
      zthres = zthres.rng[j]
      cat("Metric:",metric," | thres: top",thres,"% | outlier z <",zthres,"\n")
      res = fisher.test(merged[,metric] > quantile(merged[,metric],thres),merged$z_eOutliers < zthres)
      iter = iter + 1; res.lst[[iter]] = data.frame(iter,metric,thres,zthres,or=res$estimate,l= res$conf.int[1],h=res$conf.int[2],p = res$p.value)
    }
  }
}
res = rbindlist(res.lst)
fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_genomewide_outliers.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

fisher.test((merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020 > quantile(merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020,thres)),merged$z_eOutliers < zthres)
fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres),merged$z_eOutliers < zthres)
fisher.test((merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020 > quantile(merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020,thres)) & merged$peak_overlap.fetal_brain.Excitatory_neurons.domcke_2020,merged$z_eOutliers < zthres)
fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres) & merged$peak_overlap.Cluster1.corces_2020,merged$z_eOutliers < zthres)

fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres),merged$FLARE_ab > quantile(merged$FLARE_ab,thres))

thres=0.99;zthres = -3
table(merged$FLARE_ab > quantile(merged$FLARE_ab,thres),merged$z_eOutliers < zthres)
merged.sub = subset(merged,FLARE_ab > quantile(FLARE_ab,thres) & z_eOutliers < zthres)
merged$set = merged$id %in% merged.sub$id
table(merged$FLARE_fb > quantile(merged$FLARE_fb,0.99),merged$set)
fisher.test(merged$FLARE_fb > quantile(merged$FLARE_fb,0.99),merged$set)

fisher.test(merged$FLARE_fb > quantile(merged$FLARE_fb,thres) & merged$FLARE_ab > quantile(merged$FLARE_ab,thres),merged$z_eOutliers < zthres)














#################################

##############

library(data.table)
df.mg = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/merged/rosmap.FLARE.outliers.merged.all_chr.tsv",data.table = F,stringsAsFactors = F)
mapping = fread("/oak/stanford/groups/smontgom/amarder/data/ensg_hgnc_mapping_115.txt",data.table = F,stringsAsFactors = F)

dist_end = abs(mapping$`Transcription start site (TSS)` - mapping$`Gene end (bp)`)
dist_start = abs(mapping$`Transcription start site (TSS)` - mapping$`Gene start (bp)`)
dist_min = apply(cbind(dist_end,dist_start),1,min)
mapping = mapping[order(dist_min,decreasing = FALSE),]
mapping = mapping[!duplicated(mapping$`Gene stable ID`),]

df.mg = merge(df.mg,mapping,by.x="GeneName",by.y="Gene stable ID")
df.mg = subset(df.mg,`Gene name`!="")
dist_end = abs(df.mg$pos - df.mg$`Gene end (bp)`)
dist_start = abs(df.mg$pos - df.mg$`Gene start (bp)`)
dist_min = apply(cbind(dist_end,dist_start),1,min)
df.mg$dist_min = dist_min
dist_tss = abs(df.mg$pos - df.mg$`Transcription start site (TSS)`)
df.mg$dist_tss = dist_tss
cor.test(log10(df.mg$dist_tss+1),log10(df.mg$gene_distance_1+1))

df.mg$FLARE_diff = scale(df.mg$FLARE_ab) - scale(df.mg$FLARE_heart)
df.mg.sub = subset(df.mg,dist_tss < 2000)
merged <- df.mg.sub %>%
  left_join(df[,c("id","peak_overlap.Cluster1.corces_2020","peak_overlap.fetal_brain.Excitatory_neurons.domcke_2020","abs_logfc.mean.Cluster1.corces_2020","abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020")], by = "id")

merged$cbp = merged$abs_logfc.mean.Cluster1.corces_2020
merged$cbp_peak = ifelse(merged$peak_overlap.Cluster1.corces_2020, merged$cbp,0)
merged$cbp_nopeak = ifelse(merged$peak_overlap.Cluster1.corces_2020, 0, merged$cbp)

rng = seq(0.9,0.999,by=0.001)
zthres.rng = c(-2,-2.5,-3)
iter = 0; res.lst = list()
for (metric in c("FLARE_ab","cbp","cbp_peak","cbp_nopeak")) {
  for (j in 1:length(zthres.rng)) {
    for (i in 1:length(rng)) {
      thres = rng[i]
      zthres = zthres.rng[j]
      cat("Metric:",metric," | thres: top",thres,"% | outlier z <",zthres,"\n")
      res = fisher.test(merged[,metric] > quantile(merged[,metric],thres),merged$z_eOutliers < zthres)
      iter = iter + 1; res.lst[[iter]] = data.frame(iter,metric,thres,zthres,or=res$estimate,l= res$conf.int[1],h=res$conf.int[2],p = res$p.value)
    }
  }
}
res = rbindlist(res.lst)
fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_outliers.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



# fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),df.mg.sub$z_eOutliers > 2)
# fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.9),df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.99),df.mg.sub$z_eOutliers < -3)
table(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.99),df.mg.sub$z_eOutliers < -3)

fisher.test(merged$FLARE_ab > quantile(merged$FLARE_ab,0.99),merged$z_eOutliers < -3) # replicate above
fisher.test((merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020 > quantile(merged$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020,thres)),merged$z_eOutliers < zthres)
fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres),merged$z_eOutliers < zthres)
fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres) & merged$peak_overlap.Cluster1.corces_2020,merged$z_eOutliers < zthres)
fisher.test(merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,thres) & !merged$peak_overlap.Cluster1.corces_2020,merged$z_eOutliers < zthres)
merged.sub = merged[merged$peak_overlap.Cluster1.corces_2020,]
fisher.test(merged.sub$FLARE_ab > quantile(merged$FLARE_ab,0.99),merged.sub$z_eOutliers < zthres)
fisher.test(merged.sub$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,0.99),merged.sub$z_eOutliers < zthres)
table(merged.sub$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,0.99),merged.sub$z_eOutliers < zthres)
table(merged.sub$FLARE_ab > quantile(merged$FLARE_ab,0.99),merged.sub$z_eOutliers < zthres)

nrow(subset(df.mg.sub,merged$FLARE_ab > quantile(merged$FLARE_ab,0.99) & merged$z_eOutliers < -3))
nrow(subset(df.mg.sub,merged$abs_logfc.mean.Cluster1.corces_2020 > quantile(merged$abs_logfc.mean.Cluster1.corces_2020,0.99) & merged$z_eOutliers < -3))

topk_recall <- function(scores, labels, k) {
  ord <- order(scores, decreasing=TRUE)
  sum(labels[ord][1:k]) / sum(labels)
}

library(dplyr)

tail_enrichment <- function(scores, labels, probs = seq(0.9, 0.999, by=0.01)) {
  sapply(probs, function(p) {
    thr <- quantile(scores, p)
    tab <- table(scores > thr, labels)
    fisher.test(tab)$estimate # odds ratio
  })
}

ors <- tail_enrichment(merged$FLARE_ab, merged$z_eOutliers < -3)
plot(seq(0.9, 0.999, by=0.01), ors, type="b", xlab="Quantile cutoff", ylab="Odds ratio")


library(pROC)

library(processx)

subset(df.mg.sub,df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.99) & df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),df.mg.sub$z_eOutliers < -3)


fisher.test(df.mg.sub$PHRED > quantile(df.mg.sub$PHRED,0.99,na.rm=T),df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$PHRED > quantile(df.mg.sub$PHRED,0.999,na.rm=T),df.mg.sub$z_eOutliers < -3)

df.mg.sub = subset(df.mg,dist_tss > 2000)
setorder(df.mg.sub, +dist_tss,na.last = TRUE)
df.mg.sub <- df.mg.sub[!duplicated(df.mg.sub$id),]
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.99),df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),abs(df.mg.sub$z_eOutliers) > 3 )
fisher.test(df.mg.sub$FLARE_brain > quantile(df.mg.sub$FLARE_brain,0.99),abs(df.mg.sub$z_sOutliers) > 2 & df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$FLARE_diff > quantile(df.mg.sub$FLARE_diff,0.99),df.mg.sub$z_eOutliers < -3)

df.mg.sub
fisher.test(df.mg.sub$FLARE_diff > quantile(df.mg.sub$FLARE_diff,0.99),df.mg.sub$z_eOutliers < -3)

fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.95),abs(df.mg.sub$z_eOutliers) > 2)
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),df.mg.sub$z_eOutliers < -3)


table(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.99),df.mg.sub$z_eOutliers < -3)


fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3 & df.mg$z_pOutliers < -1)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_pOutliers < -2)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_pOutliers > 2)


fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)



table(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -3)
subset(df.mg,FLARE_ab > quantile(df.mg$FLARE_ab,0.999) & df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$z_eOutliers < -3)
fisher.test(df.mg$FLARE_heart > quantile(df.mg$FLARE_heart,0.999),df.mg$z_eOutliers < -3)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_eOutliers < -2)
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$z_eOutliers < -2)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_pOutliers < -2 & df.mg$z_eOutliers < -2)
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$z_pOutliers < -1)
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$z_pOutliers < -1)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z_pOutliers < -3 & df.mg$z_eOutliers < -3)
fisher.test(df.mg$PHRED > quantile(df.mg$PHRED,0.999,na.rm=T),df.mg$z_pOutliers < -3 & df.mg$z_eOutliers < -3)

df.mg.sub = subset(df.mg,gene_distance_1 > 2000)
fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$PHRED > quantile(df.mg.sub$PHRED,0.999,na.rm=TRUE),df.mg.sub$z_eOutliers < -3)
fisher.test(df.mg.sub$PHRED > quantile(df.mg.sub$PHRED,0.999,na.rm=TRUE),df.mg.sub$z_eOutliers < -2)
fisher.test(df.mg.sub$phylop > quantile(df.mg.sub$phylop,0.999,na.rm=TRUE),df.mg.sub$z_eOutliers < -2)

fisher.test(df.mg.sub$FLARE_brain > quantile(df.mg.sub$FLARE_brain,0.999),df.mg.sub$z_eOutliers < -2)
fisher.test(df.mg.sub$PHRED > quantile(df.mg.sub$PHRED,0.999,na.rm=TRUE),df.mg.sub$z_eOutliers < -2)

fisher.test(df.mg.sub$FLARE_ab > quantile(df.mg.sub$FLARE_ab,0.999),df.mg.sub$z_pOutliers < -2)
fisher.test(df.mg.sub$FLARE_brain > quantile(df.mg.sub$FLARE_brain,0.999),df.mg.sub$z_pOutliers < -2)

fisher.test(df.mg$PHRED > quantile(df.mg$PHRED,0.999,na.rm=T),df.mg$z_pOutliers < -1)

df.mg.sub = subset(df.mg,z_eOutliers < -2 & gene_distance_1 < 10000)
fisher.test(df.mg.sub$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg.sub$z_pOutliers < -1)

cor.test(df.mg.sub$FLARE_ab,df.mg.sub$z_pOutliers)

fisher.test(df.mg$PHRED > quantile(df.mg$PHRED,0.999,na.rm=T),df.mg$z_eOutliers < -2)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999,na.rm=T),df.mg$z_eOutliers < -2)
fisher.test(df.mg$phylop > quantile(df.mg$phylop,0.999,na.rm=T),df.mg$z_eOutliers < -2)


fisher.test(df.mg$z_pOutliers < -2,df.mg$z_eOutliers < -2)
fisher.test(df.mg$z_pOutliers > 2,df.mg$z_eOutliers > 2)



fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.99),df.mg$eOutliers > quantile(df.mg$eOutliers,0.99))
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.9),df.mg$eOutliers < quantile(df.mg$eOutliers,0.1))
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.95),df.mg$eOutliers < quantile(df.mg$eOutliers,0.1))

fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.95),abs(df.mg$eOutliers) < 0.001)


fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.01)
fisher.test(df.mg$FLARE_heart > quantile(df.mg$FLARE_heart,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.01)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.01)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.05)
fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.05)
fisher.test(df.mg$FLARE_heart > quantile(df.mg$FLARE_heart,0.999),df.mg$eOutliers < 0 & abs(df.mg$eOutliers) < 0.05)

# df.mg$z = qnorm(abs(df.mg$eOutliers))
signed_p_to_z <- function(p_signed) {
  sign(p_signed) * qnorm(pmin(1, abs(p_signed)/2), lower.tail = FALSE)
}
df.mg$z = signed_p_to_z(df.mg$eOutliers)

fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z < -1)
fisher.test(df.mg$FLARE_ab > quantile(df.mg$FLARE_ab,0.999),df.mg$z < -2)

res.lst = list()
score_thres = 0.999
score_type = "FLARE_ab"
outlier_thres = -2
j = 0
for (score_type in c("FLARE_ab"))
  for (score_thres in c(0.99,0.999)) { 
    for (outlier_thres in seq(-1,-4,by=-0.25)) {
      j = j + 1; print(j)
      res = fisher.test(df.mg[,score_type] > quantile(df.mg$FLARE_ab,score_thres),df.mg$z < outlier_thres)
      res.lst[[j]] = data.frame(score_type,score_thres,outlier_thres,or=res$estimate,l=res$conf.int[1],h=res$conf.int[2],pval=res$p.value)
    }
  }
}
as.data.frame(do.call(rbind,res.lst))


cor.test(df.mg$FLARE_brain,df.mg$z)

fisher.test(df.mg$FLARE_brain > quantile(df.mg$FLARE_brain,0.95),df.mg$eOutliers < quantile(df.mg$eOutliers,0.1))

