library(data.table)
library(ggplot2)
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/qtl_ieqtl_analysis.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$int = ifelse(df$int,"Neuron-specific","Not neuron-specific")
colnames(df)[3:4] = c("Neuron","Microglia")

g = ggplot(df,aes(x=int,y=Microglia,col=int)) + geom_jitter() +
  ggpubr::theme_pubr() +
  labs(x="Fine-mapped\nGTEx brain eQTLs",y="Microglia chromBPnet scores [ |logFC| ]",col="") +
  theme(axis.text.x = element_blank()) +
  scale_color_manual(values=c("#A8A8FF","#53D6D6")) +
  theme(legend.position = 'right')
  
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/qtl_ieqtl_analysis.pdf"
pdf(f,width = 4,height=5)
print(g)
dev.off()
