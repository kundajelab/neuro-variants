library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)

tmp = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/dev_specific_chrombpnet_motif_prop.txt",data.table = F,stringsAsFactors = F)
df_long <- tmp %>%
  pivot_longer(
    cols = starts_with("tf_"),  # Select columns that start with "tf_"
    names_to = "condition",     # Name of the new column for the variable names
    values_to = "value"         # Name of the new column for the values
  )
# df_long$value = sqrt(df_long$value)
df_long <- df_long %>%
  mutate(value = ifelse(condition %in% c("tf_fetal_specific_fb", "tf_adult_specific_fb"), -value, value))
df_long$condition = factor(df_long$condition,levels=c("tf_adult_specific_ab",
                                                      "tf_fetal_specific_ab",
                                                      "tf_adult_specific_fb",
                                                      "tf_fetal_specific_fb"))
# Create the plot
# g = ggplot(df_long, aes(x = tf, y = value), fill = condition) +
g = ggplot(df_long, aes(x = tf, y = value, fill = condition)) +
  # g = ggplot(df_long, aes(x = tf, y = sign(value) * sqrt(abs(value)), fill = condition)) +
  geom_col(position = "stack", width = 0.7,col='black') + 
  scale_fill_manual(values = c("tf_fetal_specific_fb" = "#E0DF70","tf_fetal_specific_ab" = "#E0CA70", "tf_adult_specific_ab" = "#483FA3", "tf_adult_specific_fb" = "#292361"),
                    labels=c("tf_adult_specific_ab"= "Adult-specific (ab)",
                             "tf_adult_specific_fb"= "Adult-specific (fb)",
                             "tf_fetal_specific_ab" ="Fetal-specific (ab)",
                             "tf_fetal_specific_fb" ="Fetal-specific (fb)"
                             )
                    ) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = "Transcription Factor Motif",
       y = "Proportion",
       fill = "Variant set (context)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  geom_abline(slope=0,intercept = 0,lty=1,col='black') +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(labels=abs);g

f.out = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/excitatory_neurons_dev.dataset.tf.chrombpnet.pdf"
pdf(f.out,width = 9.3,height = 6)
print(g)
dev.off()
