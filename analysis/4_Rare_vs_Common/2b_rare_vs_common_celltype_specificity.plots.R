library(data.table)
f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.celltype_specificity.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df.sub = subset(df,analysis=="orig")

g=ggplot(df.sub,aes(x=context,y=beta,ymin=beta-1.96*se,ymax=beta+1.96*se,col=context)) +
  geom_pointrange() +
  coord_flip() +
  ggpubr::theme_pubr() +
  labs(y="Mean differences (rare vs common)",
       x="Context",
       title = "# Affected Cell Types") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_color_manual(values=c("Fetal brain" = "#E0CA70",
                              "Adult brain" = "#483FA3",
                              "Fetal heart" = "#B30606",
                              "Adult heart" = "#A34D3F"
                              # "Fetal heart" = "#852222"
  )) +
  guides(col='none');g

pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/rare_vs_common.cell_type_specificity.pdf"),width = 6,height=3)
print(g)
dev.off()

