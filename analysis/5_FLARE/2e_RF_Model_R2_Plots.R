library(data.table)
f="/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/old/models/model_summary.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
modelType="rf"
df.sub = df[df$modelType==modelType,][1:2,]
df.sub$variantSet="1KG EUR Rare Variants"
# df.sub2 = df.sub;df.sub2$variantSet="ASD"
f="/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/model_summary.asd.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df.sub2 = df[df$modelType==modelType,]
df.sub2$variantSet="ASD De Novo Mutations"
df2 = as.data.frame(rbind(df.sub[,-4],df.sub2))

# Create the plot
p <- ggplot(df2, aes(x = model, y = r_squared,fill=model)) +
  geom_bar(stat = "identity",col='black') +  # Bar plot
  # geom_bar(stat = "identity", fill = "skyblue",col='black') +  # Bar plot
  geom_errorbar(aes(ymin = l^2, ymax = h^2), width = 0.2) +  # Error bars
  labs(x = "Model", y = expression("Estimated R"^2), title = "Random Forest - Model Comparison",fill="Model") +
  theme_minimal() +
  facet_wrap(.~variantSet,ncol = 2) +
  scale_fill_manual(values=c(
    
    "baseline + fb peaks" = "#F7DD4A",
    "baseline + fb peaks + cbp" =  "#E0CA70"
  ),
  labels=c(
    "baseline + fb peaks" = "Fetal brain (peaks)",
    "baseline + fb peaks + cbp" =  "Fetal brain")
  ) +
  # coord_cartesian(ylim = c(0.004,0.01)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1),
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust=0.5));p  
f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/model_r2_compare.other_mod.pdf"
pdf(f.out,width = 7,height=4)
print(p)
dev.off()
