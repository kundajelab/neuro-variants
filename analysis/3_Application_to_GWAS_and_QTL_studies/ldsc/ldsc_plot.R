library(data.table)
library(ggplot2)
library(dplyr)
f="~/Downloads/CRC_labels.csv"
f="~/Downloads/celltype_groups.tsv"
labels = fread(f,data.table = F,stringsAsFactors = F)
colnames(labels)=c("Name","Compartment")
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
labels$Compartment = firstup(labels$Compartment)
labels$Compartment[is.na(labels$Compartment)] = "Unknown"

for (analysisType in c("single-annot","union")) {
  f=paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/ldsc/results/CRC/",analysisType,"/CRC.cell_type_results.txt")

  df = fread(f,data.table = F,stringsAsFactors = F)
  # df = merge(df,labels,by="Name",all.x = TRUE)
  df = merge(df,labels,by="Name")
  df = subset(df,Name!="Unknown")
  
  # Calculate FDR using p.adjust
  df <- df %>%
    mutate(FDR = p.adjust(Coefficient_P_value, method = "fdr")) %>%
    mutate(log10_FDR = -log10(FDR))
  
  # Create the plot
  g=ggplot(df, aes(x = reorder(Name,log10_FDR,function(x){-mean(x)}), y = log10_FDR, fill = Compartment)) +
    geom_bar(stat = "identity",
             # width=1,
             col='black') +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
    ggpubr::theme_pubr() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1,color="#3D3C3CFC"),
          legend.position = "right",
          plot.title = element_text(hjust=0.5),
          plot.margin = margin(10, 10, 10, 50)) +
    labs(x = "Cell type", y = expression("-log"[10]~"(FDR)"), title = "Cell-type GWAS enrichments",
         fill = "Compartment") +
    scale_fill_manual(values=c("Stromal" = "#2AAA8A",
                                 "Immune" = "#0F52BA",
                                 "Epithelial" = "#E34234",
                               "Unknown" = "black")); g
    # scale_fill_manual(values=c("#8282C2","#D98A419E","#6BD1D1","black"));g
  f.out=paste0("~/Downloads/CRC_LDSC.",analysisType,".pdf")
  pdf(f.out,width = 13,height=6)
  print(g)
  dev.off()
}
