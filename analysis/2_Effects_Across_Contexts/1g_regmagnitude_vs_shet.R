library(data.table)
library(ggplot2)

f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/genedist_peakonly.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$set = factor(df$set,levels = c("Null","Specific","Multiple","Shared"))

library(dplyr)
df$setmax = ntile(df$s_het_1, 10)

res = list(); idx = 0; rng = seq(0.8,0.999,by=0.001)
for (thres in seq(0.8,0.999,by=0.001)) {
  print(thres)
  idx = idx + 1
  res[[idx]] = fisher.test(df$max_cbp > quantile(df$max_cbp,probs = thres), df$setmax==10)
}
or=as.numeric(unlist(lapply(res,function(x) x$estimate)))
ci=(lapply(res,function(x) x$conf.int))
l=as.numeric(unlist(lapply(ci,function(x){x[1]})))
h=as.numeric(unlist(lapply(ci,function(x){x[2]})))

res.df = data.frame(rng,or,l,h)
res.df$l[res.df$l < 0.76] = 0.76
g=ggplot(res.df,aes(x=1-rng,y=or)) + 
  geom_line() + ggpubr::theme_pubr() +
  geom_ribbon(aes(ymin=l,ymax=h),alpha=0.2,fill='grey') +
  geom_hline(yintercept=1,col='red',lty='dashed') +
  ylim(0.75,1.05) + labs(x="Regulatory Magnitude (top percentile threshold)",y=expression("Odds Ratio (" * s[het] * " > 90th percentile)"),title="Regulatory Magnitude vs Gene Constraint") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_x_log10(
    labels = function(x) paste0(x * 100, "%")
  ); print(g)
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/maxcbp_vs_shet.or.pdf"),width = 7*0.8,height=5*0.8)
print(g)
dev.off()

library(dplyr)
df$setmax = ntile(df$max_cbp, 10)
tmp = as.data.frame(cbind(aggregate(s_het_1~setmax,df,mean),aggregate(s_het_1~setmax,df,function(x){sd(x)/sqrt(length(x))}))[,c(1,2,4)])
colnames(tmp)[2:3] = c("Estimate","SE")
tmp$ci.low = tmp$Estimate - 1.96 * tmp$`SE`
tmp$ci.high = tmp$Estimate + 1.96 * tmp$`SE`
l=floor(10^(min(tmp$ci.low))/1000)*1000
h=ceiling(10^(max(tmp$ci.high))/1000)*1000
rng = seq(l,h,by=floor(((h - l) / 5)/1000)*1000)
rng = seq(40000,60000,by=5000)

# Plot with rescaled y-axis
g <- ggplot(tmp, aes(
  y = Estimate, 
  ymin = ci.low, 
  ymax = ci.high, 
  x = (10*setmax)-5,
  col = setmax
)) +
  geom_point() +
  geom_errorbar(height = 5) +
  theme_bw() +
  labs(y = expression("Mean "*italic(s)[het]), x = "Regulatory Magnitude Decile (%)",title = "Regulatory Magnitude vs Gene Constraint") +
  theme(
    # plot.title = element_text(hjust = 0.5),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 12),
    strip.background = element_rect(colour = "black", fill = "#EBF2E9"),
    strip.text = element_text(color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    # axis.text.y = element_text(hjust = 1, vjust = 1, angle = 60),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75)
  ) +
  # geom_smooth(se=F,alpha=0.9,col='lightblue') + 
  scale_color_viridis_c(option = "H") +
  scale_x_continuous(
    limits = c(0, 105), # Full range of 0-100
    breaks = seq(0, 100, 10), # Major ticks at intervals of 10
    # labels = 10*seq(0, 10, 2), # Show original setmax values as labels
    labels = c(rbind(10 * seq(0, 10, 2), rep("", length(10 * seq(0, 10, 2)) - 1)))[-12], # Show original setmax values as labels
    expand = c(0, 0) # Avoid extra padding
  ) +
  # scale_x_continuous(breaks=log10(c(rng)),labels = c(rng),limits = log10(c(min(rng)-1000,max(rng)+1000))) +
  # scale_x_continuous(breaks=log10(c(rng)),labels = c(rng),limits = log10(c(37000,57000))) +
  guides(col = "none"); g
pdf(paste0("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/FinalAnalysis/plots/maxcbp_vs_shet.pdf"),width = 7*0.8,height=4*0.8)
print(g)
dev.off()


summary(lm(max_cbp~s_het_1+gene_distance_1_log10,df))$coef[2,4]

