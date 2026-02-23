library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
df.sub = fread(f,data.table = F,stringsAsFactors = F)
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Multiple","Shared"))

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp1 >= ceiling(max(df.sub$num_peakscbp1) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp1 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp1 == 1] = "Specific"
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Multiple","Shared"))

# Fit linear model with updated num_peakscbp values for sets
summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))
res1=summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))$coef[2,]

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp05 >= ceiling(max(df.sub$num_peakscbp05) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp05 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp05 == 1] = "Specific"
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Multiple","Shared"))
res05=summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))$coef[2,]

# Annotate each row based on the number of significant peaks
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp005 >= ceiling(max(df.sub$num_peakscbp005) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp005 == 0] = "Null"
df.sub$set[df.sub$num_peakscbp005 == 1] = "Specific"
df.sub$set = factor(df.sub$set,levels=c("Null","Specific","Multiple","Shared"))
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Multiple","Shared"))
res005=summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))$coef[2,]

# Reset the "set" variable to the original:
df.sub$set = "Multiple"
df.sub$set[df.sub$num_peakscbp >= ceiling(max(df.sub$num_peakscbp) * 0.8)] = "Shared"
df.sub$set[df.sub$num_peakscbp == 0] = "Null"
df.sub$set[df.sub$num_peakscbp == 1] = "Specific"
df.sub$set = factor(df.sub$set, levels = c("Null", "Specific", "Multiple", "Shared"))
df.sub$set = factor(df.sub$set,levels = c("Null","Specific","Multiple","Shared"))
# Fit linear model with updated num_peakscbp values for sets
res01=summary(lm(gene_distance_1_log10 ~ set + s_het_1, df.sub))$coef[2,]

out = data.frame(rbind(res1,res05,res01,res005))
out$thres = c("P < 0.1","P < 0.05","P < 0.01","P < 0.005")
out

library(ggplot2)


# Calculate 95% CI
out$CI_lower <- out$Estimate - 1.96 * out$`Std..Error`
out$CI_upper <- out$Estimate + 1.96 * out$`Std..Error`

# Forest plot
g=ggplot(out, aes(x = thres, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Threshold", y = "Estimate (95% CI)", title = "Mean TSS distance of specific vs null SNPs") +
  ggpubr::theme_pubr() +
  theme(plot.title = element_text(hjust=0.5))
# pdf("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/cbp/analysis/plots/numcbp_genedist.density.pdf",width = 8*0.8,height=6*0.8)
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/variable_threshold_tss_distance.pdf"),width = 6,height=4)
print(g)
dev.off()



