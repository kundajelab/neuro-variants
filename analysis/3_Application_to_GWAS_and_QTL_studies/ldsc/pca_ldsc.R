# conda activate r

library(data.table)

ldscSetting = "single-annot"
# ldscSetting = "union"

df.lst = list()
i=0
trait_rng = c("ALS_Rheenen_2021",
              "Alzheimers_Bellenguez_2022",
              "Anorexia_nervosa_Watson_2019",
              "BMI_Yengo_2018",
              "Bipolar_Mullins_2021",
              "CAD_Aragam_2022",
              "CAD_Tcheandjieu_2022",
              "MajorDepression_Meng_2024",
              "Schizophrenia_PGCWave3_2022")

# for (annot in c("CRC","ameen_2022","corces_2020")) { 
for (annot in c("ameen_2022","domcke_2020","encode_2024","trevino_2021","corces_2020")) { 
  for (trait in trait_rng) {
    i = i + 1
    f=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/results/",annot,"/",ldscSetting,"/",trait,".cell_type_results.txt")
    df = fread(f,data.table = F,stringsAsFactors = F)
    df$trait = trait  
    df$Name = paste0(df$Name,".",annot)
    df.lst[[i]] = df
    print(trait)
    print(nrow(df))
  }
}

df=as.data.frame(do.call(rbind,df.lst))
for (trait in trait_rng) {
  ind = df$trait==trait
  df$FDR[ind] = p.adjust(df[ind,"Coefficient_P_value"],'fdr')
}
subset(df,FDR < 0.1 & trait=="Alzheimers_Bellenguez_2022")
table(df$FDR < 0.1,df$trait)


library(tidyr)

# Load necessary library
library(tidyr)
library(dplyr)

# Transform the Coefficient_P_value
df <- df %>%
  mutate(log_p_value = -log10(Coefficient_P_value))

# Reshape the data
wide_df <- df %>%
  select(Name, trait, log_p_value) %>%
  pivot_wider(names_from = Name, values_from = log_p_value) %>%
  as.data.frame()
rownames(wide_df)=wide_df[,1]; wide_df = wide_df[,-1]

# View the resulting data frame
pca = as.data.frame(prcomp(wide_df)$x)
pca$trait = rownames(pca)

library(ggplot2)
ggplot(pca,aes(x=PC1,y=PC2)) + geom_point() + geom_label(aes(label=trait))

ggplot(pca,aes(x=PC1,y=PC3)) + geom_point() + geom_label(aes(label=trait))


plot(df$Coefficient,-log10(df$Coefficient_P_value))




