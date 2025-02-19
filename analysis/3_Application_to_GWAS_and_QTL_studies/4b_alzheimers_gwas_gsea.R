library(data.table)
library(ggplot2)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/Alzheimers_GSEA.txt",data.table = F,stringsAsFactors = F)
df = subset(df,GeneCount.x >= 2 | GeneCount.y >= 2)
# cor.test(-log10(df$Adjusted.P.value.x),-log10(df$Adjusted.P.value.y),method = 'spearman')
cor.test(-log10(df$Adjusted.P.value.x),-log10(df$Adjusted.P.value.y),method = 'pearson')

g = ggplot(df,aes(x=-log10(df$Adjusted.P.value.x),y=-log10(df$Adjusted.P.value.y))) + geom_point() +
  ggpubr::theme_pubr() + 
  labs(
    x = expression("Bulk GSEA [-log"["10"]* "(" * italic(P) * ")]"),
    y = expression("Microglia-specific GSEA [-log"["10"]* "(" * italic(P) * ")]")
  );g
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/Alzheimers_GSEA.pdf"
pdf(f,width = 3.5,height = 3.5)
print(g)
dev.off()

# df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/Alzheimers_GSEA.txt",data.table = F,stringsAsFactors = F)
# df = subset(df,GeneCount.x >= 2 | GeneCount.y >= 2)
# df = df[order(df$Adjusted.P.value.y),][1:30,]
# g1 <- ggplot(df,aes(x=reorder(Term,-log10(Adjusted.P.value.y)),y=-log10(Adjusted.P.value.y),fill=Adjusted.P.value.x < 0.1)) + geom_bar(stat='identity') +
#   theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle=60,hjust=1)) +
#   labs(x='GO Term: Biological Process',y='Significance (-log10 FDR)') + scale_fill_manual(values=c('steelblue4','lightblue')) + geom_hline(yintercept = 1,col='red',lty='dashed'); g1
# g1

df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/Alzheimers_GSEA.txt",data.table = F,stringsAsFactors = F)
df = subset(df,GeneCount.x >= 2 | GeneCount.y >= 2)
df = subset(df,Adjusted.P.value.y < 0.1 & GeneCount.y >= 2)
g <- ggplot(df,aes(x=Adjusted.P.value.x > 0.1,y=-log10(Adjusted.P.value.y),col=Adjusted.P.value.x > 0.1)) +
  geom_jitter(width=0.1) + 
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,color='black'), #,colour = axis_colors),
        axis.text=element_text(color='black'),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none'
  ) +
  labs(x='GSEA Analysis',y= expression(-log[10](italic(FDR)))) + scale_color_brewer(palette="Set2") +
  scale_x_discrete(labels=c("Shared","Microglia-specific")); g

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/Alzheimers_GSEA_2.pdf"
pdf(f,width = 2.5,height = 4)
print(g)
dev.off()


df.sub = subset(df,Adjusted.P.value.x > 0.1 & Adjusted.P.value.y < 0.1)
subset(df,Adjusted.P.value.x < 0.1 & Adjusted.P.value.y < 0.1)

df.sub[order(df.sub$GeneCount.y,-1*df.sub$Adjusted.P.value.y,decreasing = T),][1:10,]

