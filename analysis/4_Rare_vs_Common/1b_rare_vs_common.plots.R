library(data.table)
library(ggplot2)
library(dplyr)

prefix="rare_vs_common"
# prefix="rare_vs_common.peaks"

for (prefix in c("rare_vs_common","rare_vs_common.peaks")) { 
  f=paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/",prefix,".txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  
  model_meta = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/model_performance.tsv",data.table = F,stringsAsFactors = F)
  model_meta <- model_meta %>%
    group_by(celltype, dataset) %>%
    mutate(
      mean_peaks_pearsonr = mean(peaks_pearsonr, na.rm = TRUE),
      mean_peaks_spearmanr = mean(peaks_spearmanr, na.rm = TRUE),
      mean_peaks_median_jsd = mean(peaks_median_jsd, na.rm = TRUE)
    ) %>%
    ungroup() %>% as.data.frame()
  model_meta = model_meta[!duplicated(paste(model_meta$celltype,model_meta$dataset)),]
  
  df = merge(df,
             model_meta[,c("dataset","celltype","dev","organ","mean_peaks_pearsonr","mean_peaks_spearmanr")],
             by.x=c("dataset","cell"),
             by.y = c("dataset","celltype"))
  df$condition = paste(df$dev,df$organ)
  
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
  
  immune_list <- c(
    "LV.lymphocyte.H", "LV.lymphocyte.I", "LV.lymphocyte.In", "LV.lymphocyte.NI",
    "fetal_heart.Myeloid_cells",
    "LV.macrophage.H", "LV.macrophage.I", "LV.macrophage.In", "LV.macrophage.NI",
    "LV.mast_cell.H", "LV.mast_cell.I", "LV.mast_cell.In", "LV.mast_cell.NI",
    "Microglia"
  )
  df$immune <- factor(ifelse(df$cell_v2 %in% immune_list, "Immune", "Non-Immune"),levels=c("Non-Immune","Immune"))
  
  # comparing fetal brain neurons to other contexts
  summary(lm(Estimate~context_v2+mean_peaks_spearmanr + immune,df))
  
  max(df$Estimate)/mean(subset(df,!fetal_neuron)$Estimate)
  df[which.max(df$Estimate),]
  
  # supp fig:
  library(ggforce)
  pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/",prefix,".pdf"),width = 8,height=15)
  g = ggplot(df, aes(x = Estimate, 
                     xmin = Estimate - 1.96 * `Std. Error`, 
                     xmax = Estimate + 1.96 * `Std. Error`, 
                     y = cell_v2,
                     col=context_v2)) +
    geom_point(aes(shape=dev)) +
    ggforce::facet_col(vars(dataset), scales = 'free_y', space = 'free') +
    # facet_grid(~dataset,scales = "free_y") +
    # facet_wrap(~context,ncol=1,scales = "free_y") +
    geom_errorbarh(height = 0.2,aes(lty=immune)) +
    theme_bw() +
    labs(x ="Diff. of pred. regulatory effects\n(rare vs common)", y = NULL,
         col="",shape="") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(colour="black",
                                        fill="#EBF2E9"),
          strip.text = element_text(color='black',face = 'bold')#,
          # panel.grid = element_blank()
    )  +
    # guides(col="none") +
    scale_shape_manual(values=c(16,17)) +
    # scale_color_viridis_d(option = "E")  +
    geom_vline(xintercept = 0,lty='dashed',col='red') +
    scale_color_manual(values=c("Fetal brain neurons" = "#E0CA70",
                                "Adult brain" = "#483FA3",
                                "Fetal brain non-neurons" = "#4D3B3B",
                                "Fetal heart" = "#FF8C69",
                                "Adult heart" = "#B30606"
    ))
  print(g)
  dev.off()
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  # main fig
  g = ggplot(df, aes(x = rank(Estimate), 
                     ymin = Estimate - 1.96 * `Std. Error`,
                     ymax = Estimate + 1.96 * `Std. Error`,
                     y = Estimate,
                     col=context_v2)) +
    geom_abline(slope = 0,intercept = 0,col='red',lty = 'dashed') +
    geom_point() +
    geom_errorbar(width = 0.2) +
    ggpubr::theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())  +
    labs(x="Rank Order",y="Diff. of pred. regulatory effects\n(rare vs common)",color="Context") +
    scale_color_manual(values=c("Fetal brain neurons" = "#E0CA70",
                                "Adult brain" = "#483FA3",
                                "Fetal brain non-neurons" = "#4D3B3B",
                                "Fetal heart" = "#FF8C69",
                                "Adult heart" = "#B30606"
    ))
  pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/",prefix,".2.pdf"),width = 8.6,height=3.3)
  print(g)
  dev.off()
  
  fwrite(df[,c("cell","context_v2","Estimate","Std. Error")],
         "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/3b.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)
  
  
  #################################
  
  
  g = ggplot(df, aes(y = Estimate,
                     x=reorder(context_v2,Estimate,median),
                     fill=context_v2)) +
    geom_violin() +
    geom_point(size=rel(0.7)) +
    coord_flip() +
    ggpubr::theme_pubr() +
    labs(y="Diff. of pred. regulatory effects\n(rare vs common)", x = "Context",
         col="",shape="") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(colour="black",
                                        fill="#EBF2E9"),
          strip.text = element_text(color='black',face = 'bold'),
          panel.grid = element_blank(),
          axis.text.x = element_text(hjust=1,vjust=1,angle=60))  +
    guides(col="none",fill="none") +
    scale_shape_manual(values=c(16,17)) +
    geom_vline(xintercept = 0,lty='dashed',col='red') +
    scale_fill_manual(values=c("Fetal brain neurons" = "#E0CA70",
                               "Adult brain" = "#483FA3",
                               "Fetal brain non-neurons" = "#4D3B3B",
                               "Fetal heart" = "#FF8C69",
                               "Adult heart" = "#B30606"
    ));
  pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/",prefix,".boxplots.pdf"),width = 4.7,height=4.7)
  print(g)
  dev.off()
  
  
}

