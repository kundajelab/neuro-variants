# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r
library(data.table)
f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/","ldsc",".FLARE.2.txt")
flare = fread(f,data.table = F,stringsAsFactors = F)

# flare = df[,c("snp_id","phylop","FLARE_baseline","FLARE_fb_peaks","FLARE_fb","FLARE_brain","FLARE_all","FLARE_heart","FLARE_ab","gene_distance_1")]

flare$id <- paste0(
  sub(":.*", "", flare$snp_id),                      # everything before first ":"
  "_",
  sub("^[^:]*:([^:]*):.*", "\\1", flare$snp_id)      # number between first and second ":"
)

for (annot_name in c("FLARE_fb","FLARE_heart")) {
  
  base_vec <- flare[[annot_name]]
  base_id  <- flare$id
  
  for (setting in c("thres95","thres99","thres999")) {

    # build tmp once per setting, impute mean once
    tmp <- data.frame(id = base_id, ANNOT = base_vec, stringsAsFactors = FALSE)
    ix_na0 <- is.na(tmp$ANNOT)
    if (any(ix_na0)) tmp$ANNOT[ix_na0] <- mean(tmp$ANNOT[!ix_na0], na.rm = TRUE)
    
    # thresholds as empirical percentiles
    if (setting != "cont") {
      p <- switch(setting,
                  "thres95"  = 0.95,
                  "thres99"  = 0.99,
                  "thres999" = 0.999)
      q <- stats::quantile(tmp$ANNOT, probs = p, na.rm = TRUE, names = FALSE)
      tmp$ANNOT <- as.integer(tmp$ANNOT > q)
    } else {
      # optional scaling for interpretability
      s <- sd(tmp$ANNOT)
      if (s > 0) tmp$ANNOT <- (tmp$ANNOT - mean(tmp$ANNOT)) / s
    }
    
    for (chrnum in 1:22) {
      message("annot=", annot_name, " | setting=", setting, " | chr ", chrnum)
      
      # df1: BIM-aligned bed (produces .bim-aligned rows)
      df1 <- fread(
        paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_bed/chr", chrnum, ".bed"),
        data.table = FALSE
      )
      # df2: lifted positions that map to your id key
      df2 <- fread(
        paste0("/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/tmp_lifted/chr", chrnum, ".lifted.bed"),
        data.table = FALSE
      )
      df2$chrpos <- paste0(df2$V1, "_", df2$V3)  # must match tmp$id scheme
      # carry SNP ID column through for matching to df1
      # (assumes df2$V4 is SNP ID that matches df1$V4)
      ix_match_df2 <- match(df2$chrpos, tmp$id)
      df2$ANNOT <- tmp$ANNOT[ix_match_df2]
      
      # any no-match NAs: fill (mean for cont, 0 for binary)
      if (setting == "cont") {
        ix_na2 <- is.na(df2$ANNOT)
        if (any(ix_na2)) df2$ANNOT[ix_na2] <- mean(df2$ANNOT[!ix_na2], na.rm = TRUE)
      } else {
        df2$ANNOT[is.na(df2$ANNOT)] <- 0L
        df2$ANNOT <- as.integer(df2$ANNOT)
      }
      
      # align to BIM order by SNP ID (thin annot requires exact row order)
      m <- match(df1$V4, df2$V4)
      annot_vec <- df2$ANNOT[m]
      
      # final NA guard post-merge
      if (setting == "cont") {
        ix_na3 <- is.na(annot_vec)
        if (any(ix_na3)) annot_vec[ix_na3] <- mean(annot_vec[!ix_na3], na.rm = TRUE)
      } else {
        annot_vec[is.na(annot_vec)] <- 0L
        annot_vec <- as.integer(annot_vec)
      }
      
      stopifnot(length(annot_vec) == nrow(df1))
      
      # THIN ANNOT OUTPUT: single column, header = annot name (or include setting)
      out <- data.frame(ANNOT = annot_vec, check.names = FALSE)
      colnames(out) <- paste0(annot_name, ".", setting)
      
      f.out <- paste0(
        "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended3/",
        annot_name, ".", setting, ".chr", chrnum, ".annot.gz"
      )
      fwrite(out, f.out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  }
}

