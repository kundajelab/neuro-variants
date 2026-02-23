library(data.table)
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/top1percent_motifbreakr.txt"
results.df = fread(f,data.table = F,stringsAsFactors = F)
head(results.df)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

df$dist_set=NA
df$dist_set[df$gene_distance_1_log10<3] = "Proximal"
df$dist_set[df$gene_distance_1_log10>3.5] = "Distal"

t.test(
  subset(df,dist_set=="Distal")$phylop,
  subset(df,dist_set=="Proximal")$phylop
)

df.sub = subset(df,set%in%c("Specific","Shared") & 
                  (df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.99))
)
table(df.sub$set)

library(ggplot2)
g = ggplot(df.sub,aes(x=gene_distance_1_log10,fill=set)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x="Distance to TSS (bp)",y="Density",fill="") +
  scale_fill_manual(
    values=c("black","#FDE725FF")
  ) +
  scale_x_continuous(
    breaks=c(0,1,2,3,4,5,6),
    labels = scales::math_format(10^.x)) + 
  theme(legend.position = c(0.45,0.95),
        legend.justification = c(1, 1))
print(g)


snps.bed.file = "~/Downloads/strongest_rv.bed"
fwrite(data.frame(df.sub$chr,df.sub$pos-1,df.sub$pos,df.sub$snp_id,0,"+"),snps.bed.file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

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
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/top1percent_motifbreakr.txt"

df.sub$dist_set=NA
df.sub$dist_set[df.sub$gene_distance_1_log10<3] = "Proximal"
df.sub$dist_set[df.sub$gene_distance_1_log10>3.5] = "Distal"
table(df.sub$dist_set,df.sub$set)

res = merge(df.sub,results.df[,c("SNP_id","geneSymbol")],by.x="snp_id",by.y="SNP_id")
res$dist_set = NA
res$dist_set[res$gene_distance_1_log10<3] = "Proximal"
res$dist_set[res$gene_distance_1_log10>3.5] = "Distal"
res = subset(res,set=="Shared")
tab = table(res$dist_set,res$geneSymbol)[c("Proximal","Distal"),]
tab = t(as.matrix(tab))
tab = as.data.frame.matrix(tab)
tab$Proximal = tab$Proximal/length(unique(res$snp_id[res$dist_set=="Proximal"]))
tab$Distal = tab$Distal/length(unique(res$snp_id[res$dist_set=="Distal"]))
tab$Enrich = log2(tab$Distal / (tab$Proximal + 2e-100))
# Distal cell-type-shared variants
tab1 = subset(tab,(Distal > 0.05 | Proximal > 0.05) & Enrich > 1.5)
tab1 = tab1[order(tab1$Distal,decreasing = T),][1:10,]
tab1
# Proximal cell-type-shared variants
tab2 = subset(tab,(Distal > 0.05 | Proximal > 0.05) & Enrich < -1.5)
tab2 = tab2[order(tab2$Proximal,decreasing = T),][1:10,]
tab2

x1 = subset(res,dist_set=="Distal" & geneSymbol %in% c("CTCF","CTCFL") & set=="Shared")
x1 = subset(df.sub,snp_id %in% x1$snp_id)
x2 = subset(res,dist_set=="Distal" & !(geneSymbol %in% c("CTCF","CTCFL")) & set=="Shared")
x2 = subset(df.sub,snp_id %in% x2$snp_id)

t.test(
  x1$max_cbp,
  x2$max_cbp
)

t.test(
  subset(df.sub,dist_set=="Distal"  & set=="Shared")$max_cbp,
  subset(df.sub,dist_set=="Proximal" & set=="Shared")$max_cbp
)

t.test(
  subset(df.sub,dist_set=="Distal" & set=="Shared")$phylop,
  subset(df.sub,dist_set=="Proximal" & set=="Shared")$phylop
)

library(ggplot2)

# Create a data frame with just the relevant columns
df_plot <- subset(df.sub, dist_set %in% c("Distal", "Proximal"))[, c("phylop", "dist_set")]

# Create the boxplot
ggplot(df_plot, aes(x = dist_set, y = phylop, fill = dist_set)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers to better see the boxes (optional)
  geom_jitter(width = 0.2, alpha = 0.3) +  # Optional: add points for individual data
  labs(
    title = "PhyloP Score Comparison",
    x = "Distance Set",
    y = "PhyloP Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



(mean(x1$max_cbp)-mean(x2$max_cbp))/mean(x2$max_cbp)


tmp = res[(res$geneSymbol=="CTCF" | res$geneSymbol=="CTCFL"),]
length(unique(tmp$snp_id))
length(unique(subset(tmp,dist_set=="Distal")$snp_id))
length(unique(subset(res,dist_set=="Distal")$snp_id))
length(unique(subset(tmp,dist_set=="Proximal")$snp_id))


length(unique(subset(res,dist_set=="Distal" & geneSymbol %in% c("CTCF","CTCFL"))$snp_id))
t.test(
  subset(res,dist_set=="Distal" & geneSymbol %in% c("CTCF","CTCFL"))$max_cbp,
  subset(res,dist_set=="Distal" & !(geneSymbol %in% c("CTCF","CTCFL")))$max_cbp
)
t.test(
  subset(res,dist_set=="Distal" & geneSymbol %in% c("CTCF","CTCFL"))$max_cbp,
  subset(res,dist_set=="Proximal")$max_cbp
)
t.test(
  subset(res,dist_set=="Distal" & !(geneSymbol %in% c("CTCF","CTCFL")))$max_cbp,
  subset(res,dist_set=="Proximal")$max_cbp
)

res = merge(df.sub,results.df[,c("SNP_id","geneSymbol")],by.x="snp_id",by.y="SNP_id")
res$dist_set = NA
res$dist_set[res$gene_distance_1_log10<3] = "Proximal"
res$dist_set[res$gene_distance_1_log10>3.5] = "Distal"
t.test(
  subset(res,set=="Shared" & dist_set=="Distal")$max_cbp,
  subset(res,set=="Specific" & dist_set=="Distal")$max_cbp
)
t.test(
  subset(res,set=="Shared" & dist_set=="Proximal")$max_cbp,
  subset(res,set=="Shared" & dist_set=="Distal")$max_cbp
)
t.test(
  subset(res,set=="Specific" & dist_set=="Distal" & geneSymbol %in% c("OLIG2"))$max_cbp,
  subset(res,set=="Specific" & dist_set=="Proximal")$max_cbp
)
t.test(
  subset(res,set=="Shared" & dist_set=="Distal" & geneSymbol %in% c("CTCF","CTCFL"))$phylop,
  subset(res,set=="Shared" & dist_set=="Proximal")$phylop
)
t.test(
  subset(res,set=="Specific" & dist_set=="Distal"  & geneSymbol %in% c("OLIG2"))$phylop,
  subset(res,set=="Specific" & dist_set=="Proximal")$phylop
)




# 50% of variants impact CTCF or CTCFL


table(res$dist_set,res$set)
tab = table(res$dist_set,res$geneSymbol)[c("Proximal","Distal"),]
tab = t(as.matrix(tab))
tab = as.data.frame.matrix(tab)
tab$Proximal = tab$Proximal/length(unique(res$snp_id[res$dist_set=="Proximal"]))
tab$Distal = tab$Distal/length(unique(res$snp_id[res$dist_set=="Distal"]))
tab$Enrich = log2(tab$Distal / (tab$Proximal + 2e-100))
# Distal cell-type-shared variants
tab1 = subset(tab,(Distal > 0.05 | Proximal > 0.05) & Enrich > 1.5)
tab1 = tab1[order(tab1$Distal,decreasing = T),][1:10,]
tab1
# Proximal cell-type-shared variants
tab2 = subset(tab,(Distal > 0.05 | Proximal > 0.05) & Enrich < -1.5)
tab2 = tab2[order(tab2$Proximal,decreasing = T),][1:10,]
tab2

library(dplyr)
library(tidyr)
results.df.tab.plot = rbind(tab1[1:10,],tab2[10:1,]) %>% as.data.frame()
results.df.tab.plot$gene_symbol = rownames(results.df.tab.plot)


