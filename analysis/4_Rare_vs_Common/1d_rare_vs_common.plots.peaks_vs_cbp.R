library(data.table)

df.lst = list()
prefix.lst = c("rare_vs_common","rare_vs_common.peaks")
for (k in 1:2) {
  prefix = prefix.lst[k]
  f=paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/",prefix,".txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  if (prefix=="rare_vs_common.diff.fetal_brain.Excitatory_neurons") {
    colnames(df)[colnames(df)=="diff"] = "Estimate"
    colnames(df)[colnames(df)=="se"] = "Std. Error"
  }
  
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
  
  df.lst[[k]] = df
}

df.mg = merge(df.lst[[2]],df.lst[[1]],by="cell")

g=ggplot(df.mg,aes(x=Estimate.x,y=Estimate.y,
                 col=context_v2.x)) +
  geom_point() +
  ggpubr::theme_pubr() +
  labs(x="Peaks Analysis",y="ChromBPnet Analysis",title="Diff. of pred. regulatory effects\n(rare vs common)",col='') +
  # guides(col="none") +
  theme(plot.title = element_text(hjust=0.5),legend.position = 'right') +
  scale_color_manual(values=c("Fetal brain neurons" = "#E0CA70",
                              "Adult brain" = "#483FA3",
                              "Fetal brain non-neurons" = "#4D3B3B",
                              "Fetal heart" = "#B30606",
                              "Adult heart" = "#A34D3F",
                              "Fetal heart" = "#852222"
  )) +
  geom_smooth(method='lm',aes(group=1),col='black',lty=3,se=F);g

pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/rare_vs_common.peaks_vs_cbp.pdf"),width = 6,height=4)
print(g)
dev.off()
cor.test(df.mg$Estimate.x,df.mg$Estimate.y)




