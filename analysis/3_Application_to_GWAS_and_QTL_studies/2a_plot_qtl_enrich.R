library(data.table)
library(ggplot2)

for (organ in c("heart_artery","brain")) {
  
  f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.qtl.",organ,".txt")
  df=fread(f,data.table = F,stringsAsFactors = F)
  # df=fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_heartarteryqtl_binarizedpip_enrichments.txt",data.table = F,stringsAsFactors = F)
  colnames(df)[5:7] = c("Std. Error","t value","Pr(<|t|)")
  model_meta = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/model_performance.tsv",data.table = F,stringsAsFactors = F)
  df = merge(df,
             model_meta[,c("dataset","celltype","dev","organ","peaks_pearsonr","peaks_spearmanr","peaks_median_jsd")],
             by.x=c("dataset","cell"),
             by.y = c("dataset","celltype"))
  df$condition = paste(df$dev,df$organ)
  # summary(lm(Estimate~condition+peaks_pearsonr+peaks_median_jsd,df))
  
  f="/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/adult_cluster_names.csv"
  labels=fread(f,data.table = F,stringsAsFactors = F)
  df = merge(df,labels[,c("Cluster","Cluster_Description")],by.x="cell",by.y="Cluster",all.x = TRUE)
  df$cell_v2 = NA; df$cell_v2[!is.na(df$Cluster_Description)] = df$Cluster_Description[!is.na(df$Cluster_Description)]
  f="/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/cell_cluster_annot.csv"
  labels=fread(f,data.table = F,stringsAsFactors = F)
  df = merge(df,labels[,c("Cluster ID","Name (long)")],by.x="cell",by.y="Cluster ID",all.x = TRUE)
  df$cell_v2[!is.na(df$`Name (long)`)] = df$`Name (long)`[!is.na(df$`Name (long)`)]
  df = df[,!(colnames(df) %in% c("Name (long)","Cluster_Description"))]
  df$cell_v2[is.na(df$cell_v2)]= df$cell[is.na(df$cell_v2)]
  
  df$neuron = grepl("neuron",df$cell_v2,ignore.case = TRUE)
  
  df$context = NA
  df$context[df$dataset=="domcke_2020" & df$neuron & df$condition=="fetal brain"] = "Fetal brain neurons (Domcke)"
  df$context[df$dataset=="domcke_2020" & !df$neuron & df$condition=="fetal brain"] = "Fetal brain non-neurons (Domcke)"
  df$context[df$dataset=="domcke_2020" & df$condition=="fetal heart"] = "Fetal heart (Domcke)"
  df$context[df$dataset=="corces_2020"] = "Adult brain (Corces)"
  df$context[df$dataset=="encode_2024"] = "Adult heart (ENCODE)"
  df$context[df$dataset=="ameen_2022"] = "Fetal heart (Ameen)"
  df$context[df$dataset=="trevino_2021" & df$neuron & df$condition=="fetal brain"] = "Fetal brain neurons (Trevino)"
  df$context[df$dataset=="trevino_2021" & !df$neuron & df$condition=="fetal brain"] = "Fetal brain non-neurons (Trevino)"
  
  df$context_v2 = sub("\\s*\\(.*", "", df$context)
  df$context_v2 = factor(df$context_v2,levels=c("Fetal brain neurons","Fetal brain non-neurons","Fetal heart","Adult heart","Adult brain"))
  
  df$fetal_neuron = df$context_v2=="Fetal brain neurons"
  
  # supp fig:
  # pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/ukb_fm_",organ,"qtl_binarizedpip_enrichments.pdf"),width = 12,height=4.8)
  pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/ukb_fm_",organ,"qtl_binarizedpip_enrichments.pdf"),width = 12,height=2.8)
  g = ggplot(df, aes(x=reorder(cell,-1*`Pr(<|t|)`,mean),
                     y = -log10(`Pr(<|t|)`), 
                     fill=context_v2)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("Fetal brain neurons" = "#E0CA70",
                               "Adult brain" = "#483FA3",
                               "Fetal brain non-neurons" = "#4D3B3B",
                               "Fetal heart" = "#B30606",
                               "Adult heart" = "#A34D3F",
                               "Fetal heart" = "#852222"
    )) +
    ggpubr::theme_pubr() +
    labs(x ="Cellular contexts (chromBPnet)", y = expression("Fine-mapped eQTL Enrichment [-log"["10"]*italic(P)*"]"),
         fill="") +
    geom_hline(col='red',lty='dashed',yintercept = -log10(0.05/nrow(df))) +
    theme(axis.text.x = element_blank())
  print(g)
  dev.off()
  
}
  
  
  