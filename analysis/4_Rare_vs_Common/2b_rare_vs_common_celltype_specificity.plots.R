library(data.table)
f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.celltype_specificity.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df.sub = subset(df,analysis=="orig")

g=ggplot(df.sub,aes(x=context,y=beta,ymin=beta-1.96*se,ymax=beta+1.96*se,col=context)) +
  geom_pointrange() +
  coord_flip() +
  ggpubr::theme_pubr() +
  labs(y="Mean differences (rare vs common)",
       x="Context",
       title = "# Affected Cell Types") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_color_manual(values=c("Fetal brain" = "#E0CA70",
                              "Adult brain" = "#483FA3",
                              "Fetal heart" = "#FF8C69",
                              "Adult heart" = "#B30606"
                              # "Fetal heart" = "#852222"
  )) +
  guides(col='none');g

pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/rare_vs_common.cell_type_specificity.pdf"),width = 6,height=3)
print(g)
dev.off()

fwrite(subset(df,analysis=="orig")[,colnames(df)!="analysis"],
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/3c.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)


library(data.table)
df.lst = list(); j=0
for (thres in c(0.005,0.01,0.05,0.1)) {
  
  f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/rare_vs_common.celltype_specificity.",thres,".txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  df$thres = thres
  df.sub = subset(df,analysis=="orig")
  j=j+1
  df.lst[[j]] = df.sub
}
df.sub = as.data.frame(do.call(rbind,df.lst))

library(ggplot2)
library(ggpubr)

g <- ggplot(df.sub,
            aes(x = context,
                y = beta,
                ymin = beta - 1.96 * se,
                ymax = beta + 1.96 * se,
                color = context,
                shape = as.factor(thres),
                linetype = as.factor(thres))) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_pubr() +
  labs(y = "Mean differences (rare vs common)",
       x = "Context",
       title = "# Affected Cell Types",
       color = NULL,shape=NULL,linetype=NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(
    "Fetal brain" = "#E0CA70",
    "Adult brain" = "#483FA3",
    "Fetal heart" = "#FF8C69",
    "Adult heart" = "#B30606"
  )) +
  scale_shape_manual(values = c(17, 16, 15, 18)) +   # pick distinct shapes for each threshold
  scale_linetype_manual(values = c("dashed", "solid", "dotdash", "dotted")) +
  guides(color = "none") # keep legend only for thres

g

pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/rare_vs_common.cell_type_specificity.thres.pdf"),width = 6,height=4.2)
print(g)
dev.off()
