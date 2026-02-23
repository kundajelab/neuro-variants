library(data.table)
library(ggplot2)
library(dplyr)

headdir="/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/"

variantSet="asd"
predictor_lst = c("FLARE_fb","FLARE_heart","PHRED","phylop","abs_logfc.mean.c11.trevino_2021","gene_distance_1.log10")
f_lst = paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/",variantSet,".",predictor_lst,".extrema.txt")
save_y_lst = list()
for (i in 1:length(f_lst)) {
  f = f_lst[i]
  predictor = predictor_lst[i]
  save_y = fread(f,data.table = F,stringsAsFactors = F)
  save_y$f = predictor
  save_y_lst[[i]] = save_y
}
save_y.all.full = as.data.frame(do.call(rbind,save_y_lst))

subset(save_y.all.full,n==20) %>% arrange(desc(prop))
save_y.all = subset(save_y.all.full,n >= 10)
save_y.all = subset(save_y.all.full,n >= 0)

# Create a subset for ymin and ymax
save_y_subset <- subset(save_y.all, f == predictor_lst[1])

# Ensure that the subset data is merged with the original data
save_y.all$ymin <- save_y_subset$l
save_y.all$ymax <- save_y_subset$h

# Make plot
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
  geom_vline(xintercept = log10(16),lty=3) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'right') +
  # ylim(0.4,0.8) +
  scale_color_manual(labels=c(FLARE_fb = "FLARE: fetal brain",
                              FLARE_heart = "FLARE: heart",
                              abs_logfc.mean.c11.trevino_2021 = "ChromBPnet: early RG",
                              PHRED = "CADD",
                              phylop = "PhyloP",
                              gene_distance_1.log10 = expression(log[10]*"(TSS distance)")),
                     values=c(FLARE_fb = "#E0CA70",
                              FLARE_heart = "#B30606",
                              abs_logfc.mean.c11.trevino_2021 = "#555599",
                              PHRED = "#48C270",
                              phylop = "#3CC2B2",
                              gene_distance_1.log10 = "black"));g

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/model_compare.ci.pdf"
pdf(f.out,width = 5.4*1.3,height=2.96*1.3)
print(g)
dev.off()

fwrite(save_y.all,
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/4c.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)



