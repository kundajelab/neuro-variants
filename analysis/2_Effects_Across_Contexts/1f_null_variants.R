library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

df = df[order(df$chr,df$pos),]

df.sub = subset(df,set%in%c("Specific") & 
                  (df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.99))
)
table(df.sub$set)

df.null = subset(df,set=="Null")

##########################################################################################

# GenomicRanges is optimized for “genomic nearest‐neighbor” queries, 
# but under the hood it’s just finding minimal absolute differences within each chromosome.
library(GenomicRanges)

# build GRanges objects
gr.null <- GRanges(
  seqnames = df.null$chr,
  ranges   = IRanges(start=df.null$pos, end=df.null$pos)
)
gr.sub  <- GRanges(
  seqnames = df.sub$chr,
  ranges   = IRanges(start=df.sub$pos, end=df.sub$pos)
)

# find the index of the nearest df.null row for each df.sub
hits <- nearest(gr.sub, gr.null)

# subset df.null by those indices
closest_df <- df.null[hits, ]

sum(df.sub$chr != closest_df$chr)
hist(df.sub$pos - closest_df$pos)

##########################################################################################

t.test(closest_df$phylop,df.sub$phylop,paired = T)
t.test(closest_df$phylop,df.sub$phylop,paired = T)$p.value
t.test(closest_df$pos,df.sub$pos,paired = T)
t.test(closest_df$pos,df.sub$pos,paired = T)$p.value

##########################################################################################

snps.bed.file = "~/Downloads/null_rv.bed"
fwrite(data.frame(closest_df$chr,closest_df$pos-1,closest_df$pos,closest_df$snp_id,0,"+"),snps.bed.file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38) 
snps.mb <- snps.from.file(snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed")

results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))
results.df = results@elementMetadata
# rownames(results.df) =NULL
results.df = as.data.frame(results.df)
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/top1percent_motifbreakr.null.txt"
fwrite(results.df,f,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

length(unique(results.df$SNP_id))
nrow(closest_df)
length(unique(results.df$SNP_id))/nrow(closest_df)
nrow(results.df)
nrow(results.df)/nrow(closest_df)
nrow(results.df)/length(unique(results.df$SNP_id))

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/top1percent_motifbreakr.txt"
results.df2 = fread(f,data.table = F,stringsAsFactors = F)
res = merge(df.sub,results.df2[,c("SNP_id","geneSymbol")],by.x="snp_id",by.y="SNP_id")

length(unique(res$snp_id))
nrow(df.sub)
length(unique(res$snp_id))/nrow(df.sub)
nrow(res)
nrow(res)/nrow(df.sub)
nrow(res)/length(unique(res$snp_id))

