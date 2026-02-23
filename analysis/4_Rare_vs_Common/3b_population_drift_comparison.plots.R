library(data.table)
library(ggplot2)
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/population_drift_comparison.txt"
df = fread(f,stringsAsFactors = F,sep='\t',data.table = F)
# df$analysis[df$pheno == "cbp_max_score"] = "Regulatory magnitude\n(max chromBPnet score)"
# df$analysis[df$pheno == "num_peakscbp"] = "# affected cell types"
# df$analysis = factor(df$analysis,levels = rev(c("# affected cell types","Regulatory magnitude\n(max chromBPnet score)")))
df$analysis[df$pheno == "cbp_max_score"] = "Regulatory magnitude"
df$analysis[df$pheno == "num_peakscbp"] = "Cell-type-specificity"
df$analysis = factor(df$analysis,levels = rev(c("Cell-type-specificity","Regulatory magnitude")))
df$popu = toupper(df$popu)
df$`Founder Population` = "No"
df$`Founder Population`[df$popu %in% c("ASJ","FIN")] = "Yes"

g=ggplot(df,aes(x=popu,col=`Founder Population`,y=Estimate,ymin=Estimate-1.96*`Std..Error`,ymax=Estimate+1.96*`Std..Error`)) +
  geom_pointrange() +
  facet_grid(.~analysis,scales="free") +
  geom_abline(slope = 0,intercept = 0,col='red',lty='dashed') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14,family="Helvetica"),
        strip.text = element_text(family="Helvetica",size=10),
        # strip.background = element_blank()) +  labs(x="Population",y="Estimated difference (SD) of MAF < 0.01 vs MAF > 0.01 SNPs",title = "Rare EUR SNPs -> Common Popu. SNPs") +
        strip.background = element_blank()) +  labs(x="Population",y="Effect size (MAF < 0.01 vs > 0.01)",title = "") +
  # guides(col="none") +
  theme(legend.position = 'top') +
  coord_flip() +
  scale_color_manual(values=c("maroon3","turquoise3"));g
print(g)

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/population_drift_comparison.plot.pdf"
pdf(f.out,width = 5,height=3)
print(g)
dev.off()


fwrite(df,
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/3d.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)

