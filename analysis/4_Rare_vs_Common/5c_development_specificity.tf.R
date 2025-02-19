library(data.table)
f = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
snps = df[,"snp_id"]

y = strsplit(snps,":")
pos = as.numeric(unlist(lapply(y,function(x){x[2]})))
chr = (unlist(lapply(y,function(x){x[1]})))

snps.bed.file = "~/Downloads/snp.bed"
fwrite(data.frame(chr,pos-1,pos,snps,0,"+"),snps.bed.file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

library(motifbreakR)
# library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
snps.mb <- snps.from.file(snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed")

results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))
results.df = results@elementMetadata
# rownames(results.df) =NULL
results.df = as.data.frame(results.df)

results.df$fetal_specific = results.df$SNP_id %in% df[df$fetal_specific,"snp_id"]

f.out = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.tf.txt"
fwrite(results.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

##########

library(data.table)
library(ggplot2)
library(tidyr)

f = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.tf.txt"
results.df = fread(f,data.table = F,stringsAsFactors = F)


# results.df = subset(results.df,SNP_id %in% subset(df,phylop > 1)$snp_id)
tab = as.data.frame.matrix(table(results.df$geneSymbol,results.df$fetal_specific))
tab$prop = tab[,2]/(tab[,1]+tab[,2])
tab$max = apply(tab[,1:2],1,max)
# length(unique(subset(results.df,fetal_specific)$SNP_id))
# length(unique(subset(results.df,!fetal_specific)$SNP_id))
# table(df$fetal_specific)
tab$prop_fetal = tab$`TRUE` / length(unique(subset(results.df,fetal_specific)$SNP_id))
tab$prop_adult = tab$`FALSE` / length(unique(subset(results.df,!fetal_specific)$SNP_id))
tab$prop_ratio = tab$prop_fetal / (tab$prop_adult+1e-30)

n_fetal = length(unique(subset(results.df,fetal_specific)$SNP_id))
n_adult = length(unique(subset(results.df,!fetal_specific)$SNP_id))

tmp = subset(tab,prop_fetal > 0.05 & prop_ratio > 4)
rownames(tmp)
tmp = subset(tab,prop_adult > 0.05 & prop_ratio < 1/4)
rownames(tmp)

# Convert to longer format
tmp = subset(tab,(prop_adult > 0.05 & prop_ratio < 1/4) | (prop_fetal > 0.05 & prop_ratio > 4))
tmp$TF = rownames(tmp)
df_long <- tmp[,c("TF","prop_fetal","prop_adult")] %>%
  pivot_longer(cols = c("prop_fetal", "prop_adult"), names_to = "Status", values_to = "Count")

# Create the plot
g = ggplot(df_long, aes(x = TF, y = Count, fill = Status)) +
  geom_col(position = "stack", width = 0.7) + 
  scale_fill_manual(values = c("prop_fetal" = "#E0CA70", "prop_adult" = "#483FA3"),
                    labels=c("prop_adult"= "Adult","prop_fetal" ="Fetal")) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = "Transcription Factor",
       y = "Proportion",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  theme(plot.title = element_text(hjust=0.5));g

f.out = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/plots/excitatory_neurons_dev.dataset.tf.pdf"
pdf(f.out,width = 8,height = 4)
print(g)
dev.off()

f.out = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.tf.txt"
fwrite(results.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



