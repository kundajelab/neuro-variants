library(data.table)
library(ggplot2)
library(dplyr)

# # Optional:
# f="~/Downloads/CRC_labels.csv"
# labels = fread(f,data.table = F,stringsAsFactors = F)
labels=NULL

trait="CAD_Tcheandjieu_2022"
ANNOTNAME="turner_2022"
for (trait in c("CAD_Tcheandjieu_2022","CAD_Aragam_2022")) {
  for (ANNOTNAME in c("IGVFCoronaryArteriesMultiome_2024","turner_2022")) { 
    for (analysisType in c("single-annot","union")) {
      f=paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/ldsc/results/",ANNOTNAME,"/",analysisType,"/",trait,".cell_type_results.txt")
      df=fread(f,data.table = F,stringsAsFactors = F)
      
      if (!is.null(labels)) {
        df = merge(df,labels,by="Name",all.x = TRUE)
      } else {
        df$Compartment = "Unknown"
      }
      
      # Calculate FDR using p.adjust
      df <- df %>%
        mutate(FDR = p.adjust(Coefficient_P_value, method = "fdr")) %>%
        mutate(log10_FDR = -log10(FDR))
      
      # Create the plot
      g=ggplot(df, aes(x = reorder(Name,log10_FDR,function(x){-mean(x)}), y = log10_FDR)) +
        ggpubr::theme_pubr() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1,color="#3D3C3CFC"),
              legend.position = "right",
              plot.title = element_text(hjust=0.5),
              plot.margin = margin(10, 10, 10, 50)) +
        labs(x = "Cell type", y = expression("-log"[10]~"(FDR)"), title = "Cell-type GWAS enrichments",
             fill = "Compartment")
      
      if (!is.null(labels)) {
        g = g + geom_bar(stat = "identity",
                         col='black',
                         aes(fill=Compartment)) +
          scale_fill_manual(values=c("#8282C2","#D98A419E","#6BD1D1","black"));
      } else {
        g = g + geom_bar(stat = "identity",
                         col='black',
                         aes(fill=log10_FDR)) +
          scale_fill_viridis_c(option='A') +
          guides(fill="none")
        # g = g + geom_bar(stat = "identity",
        #                  col='black',
        #                  fill='orange') +
        #   guides(fill="none")
      }
      g = g+geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red");g
      
      
      f.out=paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/ldsc/results/",ANNOTNAME,"/",analysisType,"/",trait,".cell_type_results.pdf")
      pdf(f.out,width = 13,height=6)
      print(g)
      dev.off()
    }
  }
}

  