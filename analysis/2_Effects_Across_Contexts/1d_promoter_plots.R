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

table(df.sub$epd,df.sub$set)

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

