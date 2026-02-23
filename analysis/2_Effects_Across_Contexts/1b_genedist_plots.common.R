library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.common.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

# df$gene_distance = (10^(df$gene_distance_1_log10)-1)
# Create the plot with custom colors and line types
g=ggplot(df, aes(x = gene_distance_1_log10, col = set, linetype = set)) +
  stat_density(geom="line",position="identity") +
  ggpubr::theme_pubr() +
  scale_color_manual(values = c("Null" = "red", "Specific" = "black", "Multiple" = "#79CDF7", "Shared" = "#FDE725FF")) +
  scale_linetype_manual(values = c("Null" = "dashed", "Specific" = "solid","Multiple" = "solid", "Shared" = "solid")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6), labels = scales::math_format(10^.x)) +
  labs(color = "Cell-type-specificity", linetype = "Cell-type-specificity",x="Distance to TSS (bp)",y="Density") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "solid", "solid"), 
                                                  color = c("red", "black", "#79CDF7", "#FDE725FF"),
                                                  size = 1),nrow=2));g
# pdf("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/cbp/analysis/plots/numcbp_genedist.density.pdf",width = 8*0.8,height=6*0.8)
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/numcbp_genedist.density.common.pdf"),width = 9.5*0.6,height=6*0.6)
print(g)
dev.off()

tmp = as.data.frame(cbind(aggregate(gene_distance_1_log10~set,df,mean),aggregate(gene_distance_1_log10~set,df,function(x){sd(x)/sqrt(length(x))}))[,c(1,2,4)])
colnames(tmp)[2:3] = c("Estimate","SE")
tmp$ci.low = tmp$Estimate - 1.96 * tmp$`SE`
tmp$ci.high = tmp$Estimate + 1.96 * tmp$`SE`
l=floor(10^(min(tmp$ci.low))/1000)*1000
h=ceiling(10^(max(tmp$ci.high))/1000)*1000
rng = seq(l,h,by=floor(((h - l) / 5)/1000)*1000)
rng = c(0,1000,5000,10000,25000,50000,75000)

g = ggplot(tmp, aes(x = Estimate, 
                    xmin = ci.low, 
                    xmax = ci.high, 
                    y = set,
                    col=set)) +
  geom_point() +
  geom_errorbarh(height = 0.4) +
  theme_bw() +
  labs(x = "Mean TSS distance (bp)",y="") +#, y = "Cell-type-\n-specificity") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background=element_rect(colour="black",
                                      fill="#EBF2E9"),
        strip.text = element_text(color='black',face = 'bold'),
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust=1,vjust=1,angle=60))  +
  guides(col="none") +
  # scale_color_viridis_d() + 
  scale_color_manual(values = c("Null" = "red", "Specific" = "black", "Multiple" = "#79CDF7", "Shared" = "#FDE725FF")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + 
  scale_x_continuous(breaks=log10(c(rng)),labels = c(rng),limits = log10(c(min(rng)-1000,max(rng)+1000)));g
# pdf("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/cbp/analysis/plots/numcbp_genedist.avg.pdf",width = 6*0.8,height=2.5*0.8)
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/numcbp_genedist.avg.common.pdf"),width = 3.5,height=2)
print(g)
dev.off()

library(dplyr)
df$setmax = ntile(df$max_cbp, 10)
tmp = as.data.frame(cbind(aggregate(gene_distance_1_log10~setmax,df,mean),aggregate(gene_distance_1_log10~setmax,df,function(x){sd(x)/sqrt(length(x))}))[,c(1,2,4)])
colnames(tmp)[2:3] = c("Estimate","SE")
tmp$ci.low = tmp$Estimate - 1.96 * tmp$`SE`
tmp$ci.high = tmp$Estimate + 1.96 * tmp$`SE`
l=floor(10^(min(tmp$ci.low))/1000)*1000
h=ceiling(10^(max(tmp$ci.high))/1000)*1000
rng = seq(l,h,by=floor(((h - l) / 5)/1000)*1000)
rng = seq(40000,70000,by=5000)

tmp

# Plot with rescaled y-axis
g <- ggplot(tmp, aes(
  x = Estimate, 
  xmin = ci.low, 
  xmax = ci.high, 
  y = (10*setmax)-5,
  col = setmax
)) +
  geom_point() +
  geom_errorbarh(height = 5) +
  theme_bw() +
  labs(x = "Mean TSS distance (bp)", y = "Decile (%)",title = "ChromBPnet Magnitude") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(colour = "black", fill = "#EBF2E9"),
    strip.text = element_text(color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 60),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75)
  ) +
  scale_color_viridis_c(option = "H") +
  scale_y_continuous(
    limits = c(0, 105), # Full range of 0-100
    breaks = seq(0, 100, 10), # Major ticks at intervals of 10
    # labels = 10*seq(0, 10, 2), # Show original setmax values as labels
    labels = c(rbind(10 * seq(0, 10, 2), rep("", length(10 * seq(0, 10, 2)) - 1)))[-12], # Show original setmax values as labels
    expand = c(0, 0) # Avoid extra padding
  ) +
  # scale_x_continuous(breaks=log10(c(rng)),labels = c(rng),limits = log10(c(min(rng)-1000,max(rng)+1000))) +
  scale_x_continuous(breaks=log10(c(rng)),labels = c(rng),limits = log10(c(41000,61000))) +
  guides(col = "none"); g
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/maxcbp_genedist.avg.common.pdf"),width = 6,height=2.5)
print(g)
dev.off()

g = ggplot(df, aes(x = df[,paste0("max_cbp")], fill = set)) +
  # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5) +
  geom_histogram(aes(y = after_stat(count / sum(count))), position = "identity", alpha = 0.5,col="black") +
  facet_wrap(~set,ncol = 1,scales = "free_y") +
  labs(x="ChromBPnet Magnitude",
       y="Frequency") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = c("Null" = "red", "Specific" = "black" , "Multiple" = "#79CDF7", "Shared" = "#FDE725FF")) +
  theme(legend.position = "none")+ # remove legend
  theme(strip.text.x = element_text(size = rel(1.1), color = "white", face = "bold"), 
        strip.background = element_rect(fill="#41416B"),
        axis.text.y = element_blank());g

pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/cbp_vs_numcbp.common.pdf"),width = 5,height=7)
print(g)
dev.off()

# filtered_data <- subset(df, df$num_peakscbp >0 & df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.999));
filtered_data <- subset(df, df$num_peakscbp >0)
g1 <- ggplot(filtered_data,aes(x=num_peakscbp,fill=set)) +
  ggpubr::theme_pubr() +
  geom_histogram(aes(y = after_stat(count / sum(count))), position = "identity", alpha = 0.5,col='black',binwidth = 1) + #,fill='#FCCC49D6') + #) +
  # xlim(-1,21) + 
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = c("Null" = "red", "Specific" = "black", "Multiple" = "#79CDF7", "Shared" = "#FDE725FF")) +
  theme(plot.title = element_text(hjust=0.5),legend.title = element_blank()) +
  labs(x="# affected cell types",y="Frequency");
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/num_cbp.common.pdf"),width = 5,height=3.3)
print(g1)
dev.off()

tab = as.data.frame(as.matrix(table(df$set)))
tab$V2 = rownames(tab)
tab = tab[,c(2,1)]
colnames(tab) = c("Set","N")
rownames(tab) = NULL
tab$Set = factor(tab$Set,levels = c("Null","Specific","Multiple","Shared"))

# Create the bar plot
g=ggplot(tab, aes(x = Set, y = log10(N), fill = Set)) +
  geom_col(colour = "black") +
  scale_y_continuous(
    breaks=c(1,2,3,4,5,6),
    labels = scales::math_format(10^.x),
    limits = c(0, max(log10(tab$N)))) + 
  # labels=c("10","100","1000","10000","100000","1000000") +
  ggpubr::theme_pubr() +
  labs(
    x = "Cell-type-specificity",
    y = "Number of Rare Variants"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Null" = "red", "Specific" = "black", "Multiple" = "#79CDF7", "Shared" = "#FDE725FF"))
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/num_cbp.set_dist.common.pdf"),width = 5,height=4)
print(g)
dev.off()


df.sub = subset(df,set%in%c("Specific","Shared") & 
                  (df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.99))
)
table(df.sub$set)

snps.bed.file = "~/Downloads/strongest_rv.bed"
fwrite(data.frame(df.sub$chr,df.sub$pos-1,df.sub$pos,df.sub$snp_id,0,"+"),snps.bed.file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

summary(lm(phylop~set,df.sub))
summary(lm(gene_distance_1_log10~set,df.sub))

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
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/specific_vs_shared.distance.common.pdf"),width = 3,height=4)
print(g)
dev.off()



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
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/top1percent_motifbreakr.common.txt"
fwrite(results.df,f,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


res = merge(df.sub,results.df[,c("SNP_id","geneSymbol")],by.x="snp_id",by.y="SNP_id")
tab = table(res$set,res$geneSymbol)[c("Specific","Shared"),]
tab = t(as.matrix(tab))
tab = as.data.frame.matrix(tab)
tab$Specific = tab$Specific/length(unique(res$snp_id[res$set=="Specific"]))
tab$Shared = tab$Shared/length(unique(res$snp_id[res$set=="Shared"]))
tab$Enrich = log2(tab$Shared / (tab$Specific + 2e-100))
tab1 = subset(tab,(Shared > 0.05 | Specific > 0.05) & Enrich > 2)
tab1 = tab1[order(tab1$Shared,decreasing = T),][1:10,]
tab2 = subset(tab,(Shared > 0.05 | Specific > 0.05) & Enrich < -2)
tab2 = tab2[order(tab2$Specific,decreasing = T),][1:10,]

library(dplyr)
library(tidyr)
results.df.tab.plot = rbind(tab1[1:10,],tab2[10:1,]) %>% as.data.frame()
results.df.tab.plot$gene_symbol = rownames(results.df.tab.plot)

# Prepare the data
results_long <- results.df.tab.plot[,c("Specific","Shared","gene_symbol")] %>%
  gather(key = "Condition", value = "Count", -gene_symbol) %>%
  mutate(Condition = factor(Condition, levels = c("Shared", "Specific")),
         Count = ifelse(Condition == "Specific", -Count, Count))

# Plot
results_long$gene_symbol = factor(results_long$gene_symbol,levels=results.df.tab.plot$gene_symbol)
g = ggplot(results_long, aes(x = gene_symbol, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "identity",width=rel(0.5),alpha=0.5) +
  geom_abline(intercept=0,slope = 0,lty='dashed',col='grey') +
  # coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(y = "Proportion", x = "Transcription Factor Motif",fill="") +
  ggpubr::theme_pubr() +
  # theme_classic() +
  # theme_minimal() +
  # guides(fill="none") +
  scale_fill_manual(values=c("Specific" = "black","Shared" = "#FDE725FF")) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1));g
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/strongest_variants.motifbreakr.common.pdf"),width = 5.6*1.2,height=3*1.2)
print(g)
dev.off()

#############

# Run EPD:

library(dplyr)
library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.common.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

df.sub = subset(df,set%in%c("Specific","Shared") & 
                  (df[,paste0("max_cbp")] > quantile(df[,paste0("max_cbp")],probs=0.99))
)

epd = fread(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/variant_input.hg38.epd.common.txt"),data.table = F,stringsAsFactors = F)
df.sub$epd = df.sub$snp_id %in% epd$V4
table(df.sub$epd,df.sub$set)
tab=table(epd=df.sub$epd,set=df.sub$set) %>% as.data.frame()
tab <- tab %>%
  group_by(set) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup() %>% as.data.frame()
# library(ggplot2)
# tab = subset(tab,set %in% c("Specific","Shared"))
# tab$set = factor(tab$set, levels = rev(levels(factor(tab$set))))
g=ggplot(subset(tab, epd == TRUE & set %in% c("Specific","Shared")), aes(x = set, y = Proportion,fill=set)) + 
  geom_bar(stat = 'identity',width = rel(0.5),col='black',alpha=0.5) + 
  ggpubr::theme_pubr() +
  # scale_x_discrete(labels=c("Shared","Specific")) +
  labs(x="Cell-type-\n-specificity",y="% of SNPs in promoters") +
  guides(fill="none") +
  scale_fill_manual(values = c("Null" = "red", "Specific" = "black", "Restrained" = "#5DC863FF","Broad" = "#79CDF7", "Shared" = "#FDE725FF")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,max(subset(tab, epd == TRUE)$Proportion)*1.05)) +
  theme(axis.title.x = element_blank())
  # +
  # coord_flip()#

# scale_x_reverse()
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/strongest_variants.epd.common.pdf"),width = 2.3,height=3.5)
print(g)
dev.off()










