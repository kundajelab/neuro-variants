library(dplyr)
library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

df.sub = subset(df,set%in%c("Specific","Shared") & 
                  (df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.99))
)

epd = fread(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/variant_input.hg38.epd.txt"),data.table = F,stringsAsFactors = F)
df.sub$epd = df.sub$snp_id %in% epd$V4

# for input into GREAT:
thres=1
thres2=0
# thres=75
# thres2=thres
options(scipen = 999)
df.sub1 = subset(df.sub,!epd & set=="Shared")
tmp = data.frame(chr=df.sub1$chr,pos0=df.sub1$pos - thres,pos1=df.sub1$pos + thres2)
fwrite(tmp,"/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/nonpromoter_shared_snv.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
df.sub2 = subset(df.sub,!epd & set=="Specific")
tmp2 <- data.frame(chr=df.sub2$chr,pos0=df.sub2$pos - thres,pos1=df.sub2$pos + thres2)
fwrite(rbind(tmp,tmp2),"/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/nonpromoter_specific_snv.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
options(scipen = 0)

# wilcox.test(df.sub1$s_het_1,df.sub2$s_het_1)
# t.test(df.sub1$s_het_1,df.sub2$s_het_1)
# 
# wilcox.test(df.sub1$phylop,df.sub2$phylop)
# t.test(df.sub1$phylop,df.sub2$phylop)
# 
# t.test(df.sub1$max_cbp,df.sub2$max_cbp)
# 
# df.sub.tmp = subset(df.sub,!epd)
# 
# summary(lm(phylop~max_cbp + set + s_het_1 + gene_distance_1_log10,df.sub.tmp))
# 
# table(df.sub$epd,df.sub$set)
# 
# run GSEA:
library(enrichR)

# Set up your gene lists
foreground_genes <- unique(df.sub1$closest_gene_1)
background_genes <- unique(df.sub2$closest_gene_1)

# Check available databases
dbs <- enrichR::listEnrichrDbs()
dbs_of_interest <- c("GO_Biological_Process_2023", "KEGG_2021_Human")

# Run enrichment for foreground genes
enrichment_results <- enrichr(foreground_genes, dbs_of_interest)

# View results for one database
head(enrichment_results[["GO_Biological_Process_2023"]])
# head(enrichment_results[["KEGG_2021_Human"]])

enrichment_results2 <- enrichr(background_genes, dbs_of_interest)
enrichment_results2 <- enrichr(background_genes, "GO_Biological_Process_2023")[[1]]
# enrichment_results2 <- enrichr(background_genes, "KEGG_2021_Human")
f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/specific_variants_enrichment.txt"
fwrite(enrichment_results2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
head(enrichment_results2[["GO_Biological_Process_2023"]],20)


tab=table(epd=df.sub$epd,set=df.sub$set) %>% as.data.frame()
tab <- tab %>%
  group_by(set) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup() %>% as.data.frame()

g=ggplot(subset(tab, epd == TRUE & set %in% c("Specific","Shared")), aes(x = set, y = 100*Proportion,fill=set)) + 
  geom_bar(stat = 'identity',width = rel(0.5),col='black',alpha=0.5) + 
  ggpubr::theme_pubr() +
  labs(x="Cell-type-\n-specificity",y="% of SNPs in promoters") +
  guides(fill="none") +
  scale_fill_manual(values = c("Null" = "red", "Specific" = "black", "Restrained" = "#5DC863FF","Broad" = "#79CDF7", "Shared" = "#FDE725FF")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,max(subset(tab, epd == TRUE)$Proportion)*1.05)) +
  theme(axis.title.x = element_blank())

# scale_x_reverse()
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/strongest_variants.epd.pdf"),width = 2.3,height=3.5)
print(g)
dev.off()

