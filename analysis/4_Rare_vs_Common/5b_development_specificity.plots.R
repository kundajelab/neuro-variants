library(data.table)
library(ggplot2)
library(cowplot)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.txt",data.table = F,stringsAsFactors = F)
head(df)
df$dev_spec = ifelse(df$adult_specific,"Adult-specific","Fetal-specific")

g1=ggplot(df,aes(x=phylop,fill=dev_spec)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x="PhyloP",y="Density",fill="") +
  scale_fill_manual(
    values=c("#483FA3","#E0CA70")
  ) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = c(1, 1))

g2=ggplot(df,aes(x=abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020,fill=dev_spec)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x="Fetal excitatory neurons chromBPnet",y="Density",fill="") +
  scale_fill_manual(
    values=c("#483FA3","#E0CA70")
  ) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = c(1, 1))

g3=ggplot(df,aes(x=abs_logfc.mean.Cluster1.corces_2020,fill=dev_spec)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x="Adult excitatory neurons chromBPnet",y="Density",fill="") +
  scale_fill_manual(
    values=c("#483FA3","#E0CA70")
  ) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = c(1, 1))

# Align plots
aligned_plots <- align_plots(g2,g3,g1, align = "v", axis = "lr")

# Combine aligned plots in a grid
p <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], aligned_plots[[3]], ncol = 1)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/excitatory_neurons_dev.dataset.plot.pdf",width = 4,height = 8)
print(p)
dev.off()

g4=ggplot(df,aes(x=s_het_1,fill=dev_spec)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x=expression("Gene Constraint (" * italic(s)[het]*")"),y="Density",fill="") +
  scale_fill_manual(
    values=c("#483FA3","#E0CA70")
  ) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = c(1, 1))

library(cowplot)

# Align plots
aligned_plots <- align_plots(g2,g3,g1,g4, align = "v", axis = "lr")

# Combine aligned plots in a grid
plot_grid(aligned_plots[[1]], aligned_plots[[2]], aligned_plots[[3]], aligned_plots[[4]], ncol = 1)

suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))

genes = unique(df$closest_gene_1[df$adult_specific])
go_result <- enrichGO(gene = genes, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",  # Specify gene identifier type
                      ont = "BP",  # Options: "BP" (Biological Process), "MF", "CC"
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2)
go_result = as.data.frame(go_result)

genes2 = unique(df$closest_gene_1[df$fetal_specific])
go_result2 <- enrichGO(gene = genes2, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",  # Specify gene identifier type
                      ont = "BP",  # Options: "BP" (Biological Process), "MF", "CC"
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,
                      universe = genes)
go_result2 = as.data.frame(go_result2)
head(go_result2)



genes = unique(df$closest_gene_1[df$adult_specific])
go_result <- enrichGO(gene = genes, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",  # Specify gene identifier type
                      ont = "BP",  # Options: "BP" (Biological Process), "MF", "CC"
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,universe = genes2)
go_result = as.data.frame(go_result)

go_result



library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.accessible.txt",data.table = F,stringsAsFactors = F)
head(df)
df$dev_spec = "Non-specific"
df$dev_spec[df$adult_specific] = "Adult-specific"
df$dev_spec[df$fetal_specific] = "Fetal-specific"

summary(lm(scale(phylop)~gene_distance_1_log10 + s_het_1 + adult_specific,df[!df$fetal_specific,]))
summary(lm(scale(phylop)~gene_distance_1_log10 + s_het_1 + fetal_specific,df[!df$adult_specific,]))
summary(lm(scale(phylop)~gene_distance_1_log10 + s_het_1 + fetal_specific,df[!df$adult_specific,]))$coef[4,4]

g1=ggplot(df,aes(x=phylop,fill=dev_spec)) + geom_density(alpha=0.3) + ggpubr::theme_pubr() +
  labs(x="PhyloP",y="Density",fill="") +
  scale_fill_manual(
    values=c("#483FA3","#E0CA70","grey")
  ) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = c(1, 1))


