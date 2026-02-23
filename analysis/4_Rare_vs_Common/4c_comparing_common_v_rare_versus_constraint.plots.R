library(data.table)
library(ggplot2)
f=paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/phylop_vs_cbp.lm.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df1 = subset(df,model==5)
f=paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/rare_vs_common.txt")
df2 = fread(f,data.table = F,stringsAsFactors = F)
df = merge(df1,df2,by=c('cell','dataset'))

model_meta = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/model_performance.tsv",data.table = F,stringsAsFactors = F)
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
           model_meta[,c("dataset","celltype","dev","organ","mean_peaks_pearsonr","mean_peaks_spearmanr","mean_peaks_median_jsd")],
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

##############################

g = ggplot(df, aes(x = Estimate.x, 
                   y = Estimate.y,
                   col=context_v2)) +
  geom_point(size=rel(1.5)) +
  # theme_bw() +
  ggpubr::theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5),
        # strip.background=element_rect(colour="black",
        #                               fill="#EBF2E9"),
        # strip.text = element_text(color='black',face = 'bold'),
        panel.grid = element_blank(),
        legend.position = 'right')  +
  labs(x="Relationship to variant constraint",y="Diff. of pred. regulatory effects\n(rare vs common)",color="Context") +
  scale_color_manual(values=c("Fetal brain neurons" = "#E0CA70",
                              "Adult brain" = "#483FA3",
                              "Fetal brain non-neurons" = "#4D3B3B",
                              "Fetal heart" = "#FF8C69",
                              "Adult heart" = "#B30606"
  ))
pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/common_vs_rare_vs_constraint.pdf",width = 6 ,height=3.5)
print(g)
dev.off()

fwrite(df[,c("Estimate.x","Estimate.y","context_v2","cell")],
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/3e.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)

cor.test(df$Estimate.x,df$Estimate.y)

###############################

# plot(df$Estimate.x,df$Estimate.y)



