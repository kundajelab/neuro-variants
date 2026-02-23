# Use 1_FLARE_Predict.R to get 'df' for ASD variants as input to this script.
# can run this script locally...

# fwrite(df,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/asd_with_flare.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

library(data.table)
library(GenomicRanges)

df = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/asd_with_flare.txt",data.table = F,stringsAsFactors = F)
# ASD mutation info
asd_info = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/Trevino_et_al_info.txt",data.table = F,stringsAsFactors = F)
asd_info$GeneSymbol = "."
asd_info$GeneSymbol[asd_info$SYMBOL!="."] = asd_info$SYMBOL[asd_info$SYMBOL!="."]
asd_info$GeneSymbol[asd_info$NEAREST!="."] = asd_info$NEAREST[asd_info$NEAREST!="."]
asd_info$Pheno = factor(asd_info$Pheno,c("control","case"))
df.mg = merge(df,asd_info[,c("snp_id","Pheno","GeneSymbol","Gene","SampleID")],by="snp_id")

# SFARI gene info
sfari = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/gene_lists/SFARI-Gene_genes_08-19-2024release_08-22-2024export.csv",data.table = F,stringsAsFactors = F)
sfari.sub = subset(sfari,syndromic==1 & (`number-of-reports` >= 10))
# sfari.sub = subset(sfari,syndromic==1 & `number-of-reports` >= 10 & `gene-score` %in% c(1,2))
# sfari gene needs to be closest TSS or closest gene body
df.mg$sfari = df.mg$closest_gene_1 %in% sfari.sub[,"gene-symbol"] | #|
  df.mg$GeneSymbol %in% sfari.sub[,"gene-symbol"]

# subset to near SFARI genes
tmp = subset(df.mg,sfari)

df.mg[1:5,c("chr","pos","snp_id")]
e2g[1:5,c('#chr',"start","end","TargetGene")]

e2g = fread("/oak/stanford/groups/smontgom/amarder/data/encode_e2g/ENCFF598EIW.bed",data.table = F,stringsAsFactors = F)


# SNPs as 1-bp ranges
gr_snp <- GRanges(
  seqnames = df.mg$chr,
  ranges   = IRanges(start = df.mg$pos, end = df.mg$pos),
  snp_id   = df.mg$snp_id
)

# e2g regions
gr_e2g <- GRanges(
  seqnames = e2g$`#chr`,
  ranges   = IRanges(start = e2g$start, end = e2g$end),
  TargetGene = e2g$TargetGene
)

# Overlap
hits <- findOverlaps(gr_snp, gr_e2g)

# Build data.frame of matches
res <- data.frame(
  snp_id     = mcols(gr_snp)$snp_id[queryHits(hits)],
  chr        = as.character(seqnames(gr_snp))[queryHits(hits)],
  pos        = start(gr_snp)[queryHits(hits)],
  TargetGene = mcols(gr_e2g)$TargetGene[subjectHits(hits)]
)

head(res)

tmp = subset(df.mg,sfari)
tmp$FLARE_fb_rk = rank(tmp$FLARE_fb)/nrow(tmp)
df2 = merge(tmp,res[,c("snp_id","TargetGene")],by="snp_id")
# df2 = merge(df.mg,res[,c("snp_id","TargetGene")],by="snp_id")
df2$genematch = df2$closest_gene_1==df2$TargetGene
df2 = df2[order(df2$genematch,decreasing = T),]
df2 = df2[!duplicated(df2$snp_id),]
table(df2$genematch)
subset(df2,FLARE_fb_rk > 0.95)[,c("FLARE_fb_rk","FLARE_fb","genematch","closest_gene_1")]

