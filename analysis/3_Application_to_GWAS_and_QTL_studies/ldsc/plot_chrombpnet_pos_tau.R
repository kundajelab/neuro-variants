library(dplyr)
library(ggplot2)

df <- data.frame(
  Disease = c("Alzheimer's", "Alzheimer's", "Schizophrenia", "Schizophrenia",
              "Alzheimer's", "Alzheimer's", "Schizophrenia", "Schizophrenia"),
  CellType = c("Microglia","Microglia","Microglia","Microglia",
               "Neurons","Neurons","Neurons","Neurons"),
  Annotation = c("Negative","Positive","Negative","Positive",
                 "Negative","Positive","Negative","Positive"),
  tau = c(4.316e-08, 6.994e-08, 2.065e-08, 5.269e-08,
          -3.527e-09, -5.361e-08, 3.281e-07, 7.607e-07),
  SE = c(9.731e-09, 1.835e-08, 2.065e-08, 5.077e-08,
         7.953e-09, 5.117e-08, 3.125e-08, 1.109e-07)
) %>%
  mutate(
    CI_low = tau - 1.96 * SE,
    CI_high = tau + 1.96 * SE
  )

print(df)

# Plot
g=ggplot(df, aes(x = Annotation, y = tau, color = CellType)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~Disease+CellType, scales = "free_y", ncol = 1) +   # force 1 column
  labs(
    x = "ChromBPNet annotation",
    y = expression(tau~"(per-SNP heritability)"),
    title = ""
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, col = 'red', lty = 'dashed') +
  scale_color_manual(values=c("black","red")) +
  guides(col="none")
print(g)
pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/ldsc/plots/chrombpnet_pos.pdf",width=3.2,height=10)
print(g)
dev.off()
