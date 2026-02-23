library(data.table)
organ="brain"
for (organ in c("brain","heart_artery")) {
  f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/qtl.",organ,".cell_type_specificity.txt")
  df.sub = fread(f,data.table = F,stringsAsFactors = F)
  df.sub$context = c("Fetal brain","Adult brain","Fetal heart","Adult heart")
  g=ggplot(df.sub,aes(x=context,y=-log10(`Pr(>|t|)`),fill=context)) +
    geom_col(col='black') +
    coord_flip() +
    ggpubr::theme_pubr() +
    labs(y="Fine-mapped Enrichment (-log10 P-values)",
         x="Context",
         title="# affected cell types vs. fine-mapping in different contexts") +
    theme(plot.title = element_text(hjust=0.5)) +
    scale_fill_manual(values=c("Fetal brain" = "#E0CA70",
                               "Adult brain" = "#483FA3",
                               "Fetal heart" = "#FF8C69",
                               "Adult heart" = "#B30606"
                               # "Fetal heart" = "#852222"
    )) +
    scale_y_continuous(expand = expansion(mult = 0,add=0.5)) +
    guides(fill='none');g
  
  pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/qtl.",organ,".cell_type_specificity.pdf"),width = 6,height=2.4)
  print(g)
  dev.off()
}
