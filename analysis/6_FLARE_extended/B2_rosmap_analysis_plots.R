library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/underexpression_outliers.txt",data.table = F,stringsAsFactors = F)
head(df)
subset(df,zthres==-3 & thres == 0.99)
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 5.5] <- 5.5

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
               col = as.factor(zthres),
               alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%"))

fwrite(df,
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/5c_1.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/underexpression_outliers.pdf",width = 4,height=4)
print(g)
dev.off()


df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/overexpression_outliers.txt",data.table = F,stringsAsFactors = F)
# df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/underexpression_outliers.txt",data.table = F,stringsAsFactors = F)
head(df)

# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 5.5] <- 5.5

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

df$zthres = factor(df$zthres,levels=rev(unique(df$zthres)))
g=ggplot(df, aes(x = (1 - thres), y = or,
               col = as.factor(zthres),
               alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("blue","red","black")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%"))

fwrite(df,
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/5c_2.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)



pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/overexpression_outliers.pdf",width = 4,height=4)
print(g)
dev.off()

subset(df,zthres==-3 & thres %in% c(0.95,0.99,0.999))
subset(df,zthres==-2 & thres %in% c(0.95,0.99,0.999))

df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_genomewide_outliers.txt",
           data.table = F,stringsAsFactors = F); thres_val = 5.5
head(df)

# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > thres_val] <- thres_val

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

df$zthres = factor(df$zthres,levels=rev(unique(df$zthres)))
g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  # facet_wrap(~ metric, scales = "free_y") +
  facet_wrap(~ metric,ncol = 1) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%"))
g

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/chrombpnet_genomewide_outliers.pdf",width = 4,height=8)
print(g)
dev.off()

df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_outliers.txt",
           data.table = F,stringsAsFactors = F); thres_val = 60

# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > thres_val] <- thres_val

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

df$zthres = factor(df$zthres,levels=rev(unique(df$zthres)))
g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  # facet_wrap(~ metric, scales = "free_y") +
  facet_wrap(~ metric,ncol = 1) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%"))
g

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/chrombpnet_promoter_outliers.pdf",width = 4,height=8)
print(g)
dev.off()

library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_genomewide_outliers.txt",data.table = F,stringsAsFactors = F)
df = subset(df,metric=="cbp")
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 5.5] <- 5.5

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "ChromBPNet Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%")) +
  ylim(0,5.5)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/chrombpnet_genomewide_outliers.cbp_only.pdf",width = 4,height=4)
print(g)
dev.off()

library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_genomewide_outliers.txt",data.table = F,stringsAsFactors = F)
df = subset(df,metric=="cbp_peak")
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 5.5] <- 5.5

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "ChromBPNet Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%")) +
  ylim(0,5.5)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/chrombpnet_genomewide_outliers.cbp_peak.pdf",width = 4,height=4)
print(g)
dev.off()


library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_outliers.txt",data.table = F,stringsAsFactors = F)
df = subset(df,metric=="FLARE_ab")
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 60] <- 60

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%")) +
  ylim(0,60)

fwrite(df,
       "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart REVISION3/SourceData/5d.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/flare_promoter_outliers.pdf",width = 4,height=4)
print(g)
dev.off()

library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_outliers.txt",data.table = F,stringsAsFactors = F)
df = subset(df,metric=="cbp")
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 60] <- 60

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%")) +
  ylim(0,60)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/cbp_promoter_outliers.pdf",width = 4,height=4)
print(g)
dev.off()


library(data.table)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis/chrombpnet_outliers.txt",data.table = F,stringsAsFactors = F)
df = subset(df,metric=="cbp_peak")
# 1) choose the zthres with largest absolute value
z_thres_of_interest <- df$zthres[which.max(abs(df$zthres))]

# (optional) cap extreme CIs
df$h[df$h > 60] <- 60

# 2) build a named alpha map for all factor levels
z_levels <- levels(factor(df$zthres))
alpha_map <- setNames(ifelse(z_levels == as.character(z_thres_of_interest), 1, 0.3),
                      z_levels)

g=ggplot(df, aes(x = (1 - thres), y = or,
                 col = as.factor(zthres),
                 alpha = as.factor(zthres))) +
  # ribbon only for the chosen |z| line
  geom_ribbon(
    data = subset(df, zthres == z_thres_of_interest),
    aes(x = (1 - thres), ymin = l, ymax = h),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.19
  ) +
  geom_line() +
  ggpubr::theme_pubr() +
  labs(
    x = "FLARE Percentile Threshold",
    y = "Odds Ratio",
    col = expression('Outlier Threshold ('*italic(z)*')')
  ) +
  # use any palette you like; viridis handles arbitrary counts nicely
  scale_color_manual(values=c("black","red","blue")) +
  # 3) apply the named alpha map so only the max-|z| line is fully opaque
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_abline(intercept = 1, slope = 0, col = "red", lty = "dashed") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10), labels = c("1%", "5%", "10%")) +
  ylim(0,60)

pdf("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/rosmap/analysis_plots/cbp_peak_promoter_outliers.pdf",width = 4,height=4)
print(g)
dev.off()


