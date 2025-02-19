library(data.table)
library(ggplot2)
library(stringr)

df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/flare_mean.txt",data.table = F,stringsAsFactors = F)
df$variantSet <- str_to_title(df$variantSet)
df$variantSet[df$variantSet=="Asd"] = "ASD"
df = subset(df,metric!="phylop")
df = subset(df,metric%in% c("FLARE_baseline","FLARE_heart","FLARE_ab","FLARE_fb"))
df$metric = factor(df$metric,levels= c("FLARE_baseline","FLARE_heart","FLARE_ab","FLARE_fb"))
df$variantSet = factor(df$variantSet,levels = c("ASD","Rare","Common"))
custom_labels=c(FLARE_baseline = "Baseline",
                FLARE_ab = "Adult brain",
                FLARE_fb =  "Fetal brain",
                FLARE_heart = "Heart")

g=ggplot(df,aes(x=variantSet,y=mu,ymin=mu-1.96*se,ymax=mu+1.96*se,col=metric)) +
  facet_wrap(metric~.,nrow = 1,labeller = labeller(metric = custom_labels)) + geom_pointrange() +
  geom_line(aes(group=metric)) +
  theme_minimal() +
  labs(y="Mean FLARE scores",x="Variant Set") +
  scale_color_manual(values=c("FLARE_baseline" = "#66BBBB",
                              "FLARE_ab" = "#483FA3",
                              "FLARE_brain" = "#BA9904",
                              "FLARE_fb_peaks" = "#F7DD4A",
                              "FLARE_fb" =  "#E0CA70",
                              "FLARE_heart" = "#B30606",
                              "FLARE_all" = "#F5B949")) +
  guides(col='none')

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/plots/flare_mean_compare.pdf"
pdf(f.out,width = 7.4,height=2.8)
print(g)
dev.off()


