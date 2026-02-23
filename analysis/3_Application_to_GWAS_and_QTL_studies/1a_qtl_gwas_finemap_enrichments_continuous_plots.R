library(data.table)
library(ggplot2)

k = 0
g = list()
for (j in 1:3) {
  if (j==1) {
    df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.qtl.best_celltype_continuous_predictions.txt",data.table = F,stringsAsFactors = F)
  } else if (j==2) {
    df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.meta_best_celltype_continuous_predictions.txt",data.table = F,stringsAsFactors = F)
  } else if (j==3) {
    df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.best_celltype_continuous_predictions.txt",data.table = F,stringsAsFactors = F)
  }
  trait_uniq = unique(df$trait)
  for (i in 1:length(trait_uniq)) {
    traitName = trait_uniq[i]
    df.sub = subset(df,trait==traitName)
    
    # Add 95% confidence intervals
    df.sub$lower <- df.sub$Estimate - 1.96 * df.sub$se
    df.sub$upper <- df.sub$Estimate + 1.96 * df.sub$se
    
    # Plot
    k = k + 1
    g[[k]] = ggplot(df.sub, aes(x = thres, y = Estimate)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      labs(
        x = "PIP Threshold",
        y = "ChromBPNet Diff.",
        title = traitName
      ) +
      xlim(0.09,0.91) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust=0.5))
    
    f.out = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/continuous.",traitName,".pdf")
    pdf(f.out,width = 5,height = 3)
    print(g[[k]])
    dev.off()
  }
}

g[[1]]
g[[2]]
library(cowplot)
plot_grid(g)

library(patchwork)

# Display all plots in a grid
f.out = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/continuous.all.pdf")
pdf(f.out,width = 6,height = 9)
wrap_plots(g,ncol = 2)
dev.off()

