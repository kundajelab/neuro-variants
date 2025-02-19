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

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/results/SUMMARIZED_DATA/SUMMARIZED_DATA.txt"
fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# subset(df,FDR < 0.1 & trait=="Alzheimers_Bellenguez_2022")
subset(df,FDR < 0.1 & trait=="ALS_Rheenen_2021")
table(df$FDR < 0.1,df$trait)


# Load necessary library
library(data.table)
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/ldsc/results/SUMMARIZED_DATA/SUMMARIZED_DATA.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
df$log10_FDR = -log10(df$FDR)
# mtx = df[,c("trait","Name","FDR")]
df$sig = "Null"
df$sig[df$Coefficient_P_value < 0.05] = "Nominal"
df$sig[df$FDR < 0.1] = "Sig"
df$sig[df$FDR < 0.05] = "Sig05"
df$sig[df$FDR < 0.01] = "Sig01"
df$sig = factor(df$sig,levels=c("Null","Nominal","Sig","Sig05","Sig01"))

############

meta = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/SuppTab/SuppTable1.txt",data.table = F,stringsAsFactors = F)
df$dataset <- sub("\\..*", "", df$Name)
df$cell <- sub("^[^.]+\\.(.+)\\.[^.]+$", "\\1", df$Name)
meta = merge(meta,df,by=c("dataset","cell"))

############

meta$context = paste(meta$dev,meta$cell)
meta$cell_v3 = sub("fetal_brain\\.|fetal_heart\\.","",meta$cell_v2)
meta = meta[order(meta$dataset,meta$context,meta$cell_v3),]
meta$cell_v3 = factor(meta$cell_v3,unique(meta$cell_v3))
meta$context = paste(meta$dev,meta$organ)
g=ggplot(meta, aes(x = trait, y = cell_v3, fill = sig)) +
  geom_tile(color = "white") + # Add tile borders
  scale_fill_manual(values=#c("Null" = "#FF5B4C",
                      #   "Nominal" = "#045685",
                      # "Sig" = "#1379BD",
                      # "Sig05" = "#1B91E0",
                      # "Sig01" = "#87CFFFFF"),
                      c("Null" = "#87CFFFFF",
                        "Nominal" = "#C75D20",
                        "Sig" = "#FFC04C",
                        "Sig05" = "#FF8A7FFF",
                        "Sig01" = "#F7210A"),
                    labels=c("Null" = "Not Sig",
                             "Nominal" = "P < 0.05",
                             "Sig" = "FDR < 0.1",
                             "Sig05" = "FDR < 0.05",
                             "Sig01" = "FDR < 0.01")) +
  # scale_fill_gradient(low = "blue", high = "red") + # Customize colors
  theme_minimal() + # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = rel(0.6)), # Rotate x-axis labels
    axis.text.y = element_text(size=rel(0.6)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(fill="Significance") +
  facet_wrap(~dataset+context, scales = "free_x",nrow = 1);g

f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/ldsc/results/SUMMARIZED_DATA/heatmap.pdf"
pdf(f.out,width = 8.7,height = 8)
print(g)
dev.off()


subset(df,FDR < 0.1 & trait=="Anorexia_nervosa_Watson_2019")



