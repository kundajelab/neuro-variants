library(data.table)

f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/scores_cntnap2.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df = df[,-1]
wilcox.test(df$abs_logfc.mean.c11.trevino_2021~df$Pheno)

colnames(df) = c("Pheno","FLARE-fb","FLARE-h","PhyloP","CADD","Early RG cbp")

library(tidyr)
df_long <- pivot_longer(df, cols = -Pheno, names_to = "Feature", values_to = "Value") %>% as.data.frame()

# Create scatterplots with jittered points using ggplot2
g = ggplot(df_long, aes(x = Pheno, y = Value, color = Pheno)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 2) + 
  facet_wrap(~ Feature, scales = "free_y", nrow = 1) + 
  labs(title = "Non-coding variants near CNTNAP2", x = "Case Status", y = "Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = 'none') +
  scale_color_manual(values = c("red", "black"))

# Print the plot
print(g)

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/scores_cntnap2.pdf"
pdf(f.out,width = 4.81*1.3,height=2.96*1.3)
print(g)
dev.off()


f = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/scores_sfari.txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df = df[,-1]
t.test(df$abs_logfc.mean.c11.trevino_2021~df$Pheno)
t.test(df$FLARE_fb~df$Pheno)

colnames(df) = c("Pheno","FLARE-fb","FLARE-h","PhyloP","CADD","Early RG cbp")

library(tidyr)
df_long <- pivot_longer(df, cols = -Pheno, names_to = "Feature", values_to = "Value") %>% as.data.frame()

# Create scatterplots with jittered points using ggplot2
g = ggplot(df_long, aes(x = Pheno, y = Value, color = Pheno)) +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.8, size = 0.8) + 
  facet_wrap(~ Feature, scales = "free_y", nrow = 1) + 
  labs(title = "Non-coding variants near syndromic ASD genes", x = "Case Status", y = "Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = 'none') +
  scale_color_manual(values = c("red", "black"))

# Print the plot
print(g)

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/scores_sfari.pdf"
pdf(f.out,width = 4.81*1.3,height=2.96*1.3)
print(g)
dev.off()

