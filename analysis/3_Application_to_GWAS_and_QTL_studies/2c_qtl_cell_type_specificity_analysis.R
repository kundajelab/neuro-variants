# srun --account=smontgom --partition=batch --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r


library(data.table)
library(stringr)

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)
common$chr_pos = paste0(common$chr,"_",common$pos)

organ = "brain"
# organ = "heart_artery"
finemap = fread(paste0("/oak/stanford/groups/smontgom/amarder/data/finemap_data/gtex/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.",organ,".txt.gz"),data.table = F,stringsAsFactors = F)
finemap$chr_pos <- sub("(_[A-Z]+_[A-Z]+_b38)", "", finemap$V5)
finemap = finemap[order(finemap$V6,decreasing = T),]
finemap = finemap[!duplicated(finemap$chr_pos),]

common$gwas = NA
finemap.sub = subset(finemap,V6 > 0.9)
idx = common$chr_pos %in% finemap.sub$chr_pos
common$gwas[idx] = TRUE
finemap.sub = subset(finemap,V6 < 0.01)
idx = common$chr_pos %in% finemap.sub$chr_pos
common$gwas[idx] = FALSE
table(common$gwas)

df <- subset(common, gwas) # brain QTLs
# df2 <- subset(common, !gwas) # non fine mapped

common$gwas2 = FALSE; common$gwas2[common$gwas] = TRUE
thres = 0.01
tmp = subset(common,(abs_logfc.mean.pval.Cluster1.corces_2020 < thres | abs_logfc.mean.pval.Cluster24.corces_2020 < thres))
tmp$shared = tmp$abs_logfc.mean.pval.Cluster1.corces_2020 < thres & tmp$abs_logfc.mean.pval.Cluster24.corces_2020 < thres
summary(mod <- glm(shared ~ abs_logfc.mean.Cluster1.corces_2020 + abs_logfc.mean.Cluster24.corces_2020 + gnomADg_AF_joint_nfe + gwas2 + gene_distance_1 + s_het_1,tmp,family=binomial(link="logit")))
table(tmp$gwas,tmp$shared)

##########################

# ieqtl analysis:


df <- subset(common, gwas) # brain specific

brain_ieqtl = fread("/oak/stanford/groups/smontgom/amarder/data/finemap_data/gtex/GTEx_Analysis_v8_ieQTL/Brain_Cerebellum.Neurons.ieQTL.eigenMT.annotated.txt.gz",data.table = F,stringsAsFactors = F)
brain_ieqtl$chr_pos <- sub("(_[A-Z]+_[A-Z]+_b38)", "", brain_ieqtl$variant_id)

df$int = df$chr_pos %in% brain_ieqtl$chr_pos

table(df$chr_pos %in% brain_ieqtl$chr_pos)

t.test(abs_logfc.mean.Cluster1.corces_2020~int,df)
t.test(abs_logfc.mean.Cluster24.corces_2020~int,df)
# summary(lm(abs_logfc.mean.Cluster24.corces_2020~int + s_het_1 + log10(gene_distance_1 + 1),df))

f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/qtl_ieqtl_analysis.txt")
fwrite(df[,c("snp_id","int","abs_logfc.mean.Cluster1.corces_2020","abs_logfc.mean.Cluster24.corces_2020")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
