library(data.table)
library(stringr)
library(ggplot2)
df2 = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/flare_performance.txt",data.table = F,stringsAsFactors = F)
df2$variantSet <- str_to_title(df2$variantSet)
df2$variantSet[df2$variantSet=="Asd"] = "ASD"
df2$variantSet = factor(df2$variantSet,levels=c("Common","Rare","ASD"))
df2$model = factor(df2$model,levels=c("FLARE_baseline",
                                      "FLARE_ab" ,
                                      "FLARE_heart" ,
                                      "FLARE_fb_peaks" ,
                                      "FLARE_fb",
                                      "FLARE_brain",
                                      "FLARE_all"))



# Create the plot
p <- ggplot(df2, aes(x = model, y = r,fill=model)) +
  geom_bar(stat = "identity",col='black') +  # Bar plot
  # geom_bar(stat = "identity", fill = "skyblue",col='black') +  # Bar plot
  geom_errorbar(aes(ymin = l, ymax = h), width = 0.2) +  # Error bars
  # labs(x = "Model", y = expression("Estimated R"^2), title = "Model Comparison",fill="Model") +
  labs(x = "Model", y = expression("PhyloP Corr. (Pred vs Obs)"), title = "Model Comparison",fill="Model") +
  theme_minimal() +
  # facet_wrap(.~variantSet,ncol = 2,scales="free_y") +
  facet_wrap(.~variantSet,ncol = 3) +
  scale_fill_manual(values=c("FLARE_baseline" = "#66BBBB",
                             "FLARE_ab" = "#483FA3",
                             "FLARE_brain" = "#BA9904",
                             "FLARE_fb_peaks" = "#F7DD4A",
                             "FLARE_fb" =  "#E0CA70",
                             "FLARE_heart" = "#B30606",
                             "FLARE_all" = "#F5B949"
  ),
  labels=c("FLARE_baseline" = "Baseline",
           "FLARE_ab" = "Adult brain",
           "FLARE_brain" = "Adult & fetal brain",
           "FLARE_fb_peaks" = "Fetal brain (peaks)",
           "FLARE_fb" =  "Fetal brain",
           "FLARE_heart" = "Heart",
           "FLARE_all" = "All brain & heart")
  ) +
  # coord_cartesian(ylim = c(0.004,0.01)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1),
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust=0.5));p  
f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/model_r_compare.pdf"
pdf(f.out,width = 7*0.8,height=4*0.8)
print(p)
dev.off()
