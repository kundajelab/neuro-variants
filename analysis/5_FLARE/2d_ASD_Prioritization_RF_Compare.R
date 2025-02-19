library(data.table)

library(ggplot2)
save_y.all = list()
save_y.all[[1]] = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/3.baseline + fb peaks + cbp.asd.rf.extrema.txt",data.table = F,stringsAsFactors = F)
save_y.all[[1]]$f = "Fetal brain (rf)"
save_y.all[[2]] = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/asd.FLARE_fb.extrema.txt",data.table = F,stringsAsFactors = F)
save_y.all[[2]]$f = "Fetal brain (lasso)"
save_y.all.full = as.data.frame(do.call(rbind,save_y.all))

save_y.all = subset(save_y.all.full, n <= 1000)
save_y.all$f = factor(save_y.all$f,levels = unique(save_y.all$f))
# Create a subset for ymin and ymax
save_y_subset <- subset(save_y.all, f == "Fetal brain (lasso)")

# Ensure that the subset data is merged with the original data
save_y.all$ymin <- save_y_subset$l
save_y.all$ymax <- save_y_subset$h

# which.min(save_y_subset$pval)
# save_y_subset[35,]
# subset(save_y.all,n==44)
# save_y_subset[35,]
g = ggplot(save_y.all,
           aes(x=log10(n),
               y=prop)
) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = '#E3E0E0', alpha = 0.4) +
  geom_line(aes(col=f)) +
  geom_abline(slope=0,intercept = 0.4971098,col='red',lty='dashed') +
  ggpubr::theme_pubr() +
  labs(x="Outlier Variant Rank\n(near syndromic ASD genes)",y="Proportion of case mutations",title="",col=NULL) +
  scale_x_continuous(
    breaks = c(1, 2, 3),          # Positions of the breaks
    labels = c("10", "100","1000")    # Corresponding labels 
  ) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'right') +
  # ylim(0.4,0.8) +
  scale_color_manual(labels=c("Fetal brain (lasso)" ="Fetal brain (lasso)",
                              "Fetal brain (rf)" ="Fetal brain (rf)"
  ),
  values=c("Fetal brain (lasso)" ="#E0CA70",
           "Fetal brain (rf)" ="#BFAAAA"));g

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/model_compare.ci.rf.pdf"
pdf(f.out,width = 4.81*1.3,height=2.96*1.3)
print(g)
dev.off()


