library(data.table)

# ################################################################################
# 
# # model list
# model_lst = c("baseline",
#               "baseline + fb peaks",
#               "baseline + fb peaks + cbp",
#               "baseline + brain peaks + cbp",
#               "complete",
#               "baseline + heart peaks + cbp",
#               "baseline + ab peaks + cbp"
# )
# variantSet="rare"
# i = 5
# # df = initial_data_load(variantSet)
# # Filter SNPs with missing PhyloP values
# cat("Filter SNPs with missing PhyloP values... \n")
# ind = !is.na(df$phylop)
# df = df[ind,]
# 
# cat("Performing misc. filtering... \n")
# 
# # Use log10 distance to nearest TSS in the model:
# df$gene_distance_1.log10 = log10(df$gene_distance_1+1)
# 
# # Create dummy variables for the categorical variable
# ind = !grepl("splice",df$Consequence)
# cat("Removing ",nrow(df) - sum(ind)," splice-related SNPs...\n")
# df = df[ind,]
# 
# out = make_FLARE_input_data(variantSet,i,df)
# df2 = out[[1]]
# cols = out[[2]]
# rm(out)
# length(which(cols %in% colnames(df2)))
# x = as.matrix(df2[,cols])
# 
# # Compute standard deviations of features
# feature_sds <- apply(x, 2, sd)
# 

# Extract raw coefficients
model_lst = c("baseline",
              "baseline + fb peaks",
              "baseline + fb peaks + cbp",
              "baseline + brain peaks + cbp",
              "complete",
              "baseline + heart peaks + cbp",
              "baseline + ab peaks + cbp"
)
modelNum=3
final_mod.all.lst = list()
for (modelNum in 1:length(model_lst)) {
  print(modelNum)
  final_mod = list()
  model = model_lst[modelNum]
  for (chrNum in 1:22) {
    # Load model
    f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/models/",modelNum,".",model,".chr",chrNum,".lasso.rds")
    final_mod[[chrNum]] = readRDS(f)
  }
  
  final_mod.all = do.call(cbind,lapply(final_mod,function(x) x$beta))
  # final_mod.all = data.frame(feature=rownames(final_mod.all),weight=log10(1+apply(final_mod.all,1,mean)))
  final_mod.all = data.frame(feature=rownames(final_mod.all),weight=(apply(final_mod.all,1,mean)))
  rownames(final_mod.all) = NULL
  
  library(dplyr)
  library(stringr)
  # Filter only features with peak overlaps (or whatever pattern you're interested in)
  peak_features <- final_mod.all %>%
    filter(grepl("^peak_overlap\\.", feature))  %>%
    mutate(cell_type = str_extract(feature, "(?<=peak_overlap\\.).*"))
  peak_features= peak_features$cell_type
  final_mod.all[grep("^int_",final_mod.all[,1]),1] = paste0("int.",peak_features)
  final_mod.all.lst[[modelNum]] = final_mod.all
}
final_mod.all.lst[[1]]

tmp = merge(final_mod.all.lst[[5]],final_mod.all.lst[[6]],by='feature',all=T)



# Add identifier to each data frame before merging
final_mod.all.lst <- lapply(seq_along(final_mod.all.lst), function(i) {
  df <- final_mod.all.lst[[i]]
  colnames(df)[2] <- paste0("weight_", model_lst[i])
  df
})

# Merge all by "feature"
merged_df <- Reduce(function(x, y) full_join(x, y, by = "feature"), final_mod.all.lst)
fwrite(merged_df,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/flare/flare_lasso_weights.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

###############################################################

library(data.table)
library(ggplot2)
library(dplyr)

merged_df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/flare/flare_lasso_weights.txt",data.table = F,stringsAsFactors = F)
merged_df[,2:ncol(merged_df)] = merged_df[,2:ncol(merged_df)] != 0
# merged_df[1:100,1:4]
merged_df[,1] = sub("abs_logfc.mean","abs_logfc_mean",merged_df[,1])
x = grep("peak_overlap.", merged_df[,1], value = TRUE)
x  = merged_df[,1]
# Extract dataset (everything after last period)
merged_df$dataset <- str_extract(x, "[^.]+$")  # => "domcke_2020"

# Extract cell type (everything between first and last period)
merged_df$cell <- str_match(x, "^[^.]+\\.(.*)\\.[^.]+$")[,2]

merged_df$type=NA
merged_df$type[merged_df$feature %in% c("s_het_1","gene_distance_1.log10")] = "base"
merged_df$type[grep("peak_overlap",merged_df$feature)] = "peak"
merged_df$type[grep("abs_logfc_mean",merged_df$feature)] = "cbp"
merged_df$type[grep("int\\.",merged_df$feature)] = "int"


df = merged_df

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
           by.y = c("dataset","celltype"),all.x=T)
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


> df[1,]
cell    dataset                     feature weight_baseline weight_baseline + fb peaks weight_baseline + fb peaks + cbp weight_baseline + brain peaks + cbp
1  aCM ameen_2022 peak_overlap.aCM.ameen_2022              NA                         NA                               NA                                  NA
weight_complete weight_baseline + heart peaks + cbp weight_baseline + ab peaks + cbp type   dev organ mean_peaks_pearsonr mean_peaks_spearmanr   condition
1           FALSE                               FALSE                               NA peak fetal heart           0.7759489             0.755644 fetal heart
cell_v2 neuron             context  context_v2 fetal_neuron
1     aCM  FALSE Fetal heart (Ameen) Fetal heart        FALSE

library(ggplot2)
library(reshape2)
library(dplyr)

df_long <- melt(df, 
                measure.vars = paste0("weight_", model_lst),
                variable.name = "model", 
                value.name = "status")

# Clean up model names
df_long <- df_long %>%
  mutate(model_clean = gsub("weight_", "", model),
         model_clean = factor(model_clean, levels = model_lst))

# Plot with feature as y-axis but hide the labels
# df_long$model_clean = factor(df_long$model_clean,levels = c("baseline" ,"baseline + fb peaks","baseline + fb peaks + cbp", "baseline + ab peaks + cbp","baseline + brain peaks + cbp","baseline + heart peaks + cbp","complete" ))
# Define models and their pretty labels
model_lst <- c(
  "baseline",
  "baseline + fb peaks",
  "baseline + fb peaks + cbp",
  "baseline + brain peaks + cbp",
  "complete",
  "baseline + heart peaks + cbp",
  "baseline + ab peaks + cbp"
)

pretty_labels <- c(
  "Baseline",
  "Fetal brain (peaks)",
  "Fetal brain",
  "Adult & fetal brain",
  "All brain & heart",
  "Heart",
  "Adult brain"
)

# Clean and label model column
df_long$model_clean <- factor(
  gsub("weight_", "", df_long$model),
  levels = model_lst,
  labels = pretty_labels
)
# df_long$model_clean = factor(df_long$model_clean,levels = c(
#   "Baseline",
#   "Fetal brain (peaks)",
#   "Fetal brain",
#   "Adult brain",
#   "Adult & fetal brain",
#   "Heart",
#   "All brain & heart"
# ))


df_long$type = factor(df_long$type,levels = c("base","peak","cbp","int"))
df_long$feature = factor(df_long$feature,levels=unique(df_long[order(df_long$type,df_long$context_v2),]$feature))
g=ggplot(df_long, aes(x = model_clean, y = feature)) +
  # Background grey dot for FALSE or NA
  # geom_point(data = subset(df_long, status != TRUE),
  #            color='white',fill='white', shape = 21, size = 3, stroke = 0.2) +
  # Colored dot for TRUE
  geom_point(data = subset(df_long, status == TRUE), 
             # aes(color = context_v2,fill = context_v2), shape = 21, size = 3, stroke = 0.2) +
             aes(color = context_v2,fill = context_v2,shape=type), size = 3, stroke = 0.2) +
  # facet_wrap(~type, ncol = 1,) +
  scale_color_manual(values = c(
    "Fetal brain neurons" = "#E0CA70",
    "Adult brain" = "#483FA3",
    "Fetal brain non-neurons" = "#4D3B3B",
    "Fetal heart" = "#B30606",
    "Adult heart" = "#A34D3F",
    "Fetal heart" = "#852222"
  )) + 
  scale_fill_manual(values = c(
    "Fetal brain neurons" = "#E0CA70",
    "Adult brain" = "#483FA3",
    "Fetal brain non-neurons" = "#4D3B3B",
    "Fetal heart" = "#B30606",
    "Adult heart" = "#A34D3F",
    "Fetal heart" = "#852222"
  )) +
  scale_shape_manual(
    values = c(
      "base" = 19,
      "peak" = 19,
      "cbp" = 15,
      "int" = 17
    ),
    labels=c("Base","Peaks","ChromBPNet","ChromBPNet x Peak")
  ) +
  # theme_minimal() +
  theme_bw() +
  # ggpubr::theme_pubr() +
  theme(
    axis.text.y = element_blank(),  # <- hides feature names
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust=0.5)
    # panel.grid = element_blank(),
  ) +
  labs(x = "Model", y = NULL, color = "Context",shape="Feature Set",
       title = "Non-Zero FLARE Coefficents") +
  guides(fill="none")

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/flare/flare_lasso_weights.pdf",width = 5,height=8)
print(g)
dev.off()

