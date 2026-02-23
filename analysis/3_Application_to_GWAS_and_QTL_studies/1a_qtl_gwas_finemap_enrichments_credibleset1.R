trait_lst = c("BMI","Neuroticism","AFib")
traitName="BMI"
traitName="Alzheimer_LTFH"
j = 0; res.lst = list()

f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.",traitName,".txt")
res.df = fread(f,data.table = F,stringsAsFactors = F)
res.df = res.df[which.min(res.df$pval),]

print(traitName)
finemap.ukb.sub = subset(finemap.ukb,trait==traitName)
finemap.ukb.sub$chr_pos = paste0(finemap.ukb.sub$chromosome,"_",finemap.ukb.sub$end_hg38)

finemap.ukb.sub = subset(finemap.ukb.sub,cs_id != -1)
finemap.ukb.sub$cs = apply(finemap.ukb.sub[,c("region","method","cs_id")],1,paste,collapse="-")

length(unique(finemap.ukb.sub$region))
length(unique(finemap.ukb.sub$cs))
length(unique(subset(finemap.ukb.sub,method=="FINEMAP")$cs))
length(unique(subset(finemap.ukb.sub,method=="SUSIE")$cs))

tab = table(finemap.ukb.sub$cs)
tab = tab[tab==1]
length(tab)
finemap.ukb.sub.cs1 = subset(finemap.ukb.sub,cs %in% names(tab))
length(unique(finemap.ukb.sub.cs1$chr_pos))
length(unique(finemap.ukb.sub.cs1$region))
common.sub = subset(common,chr_pos %in% finemap.ukb.sub.cs1$chr_pos)
nrow(common.sub)


string = grep("abs_logfc.mean.pval.",colnames(common),value = TRUE)
dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")
#
cell = res.df$cell[1]
dataset = res.df$dataset[1]

cell = cell_lst[31]
dataset = dataset_lst[31]
common.sub[,paste0("peak_overlap.",cell,".",dataset)]
common.sub[,paste0("abs_logfc.mean.",cell,".",dataset)]
common.sub[,paste0("abs_logfc.mean.pval.",cell,".",dataset)]


traitName = "Alzheimers_Bellenguez_2022"
print(traitName)
fm = fread(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/all/",traitName,".susie_all.txt"),data.table = F,stringsAsFactors = F)
fm.sub = subset(fm,CS.L_10 != -1 | CS.L_1 != -1)
tmp = subset(fm.sub,CS.L_10!= -1)[,1:5]
colnames(tmp)[3:5] = c("converged","CS","PIP")
tmp$method = 10
tmp1 = tmp
tmp = subset(fm.sub,CS.L_1!= -1)[,c(1:2,6:8)]
colnames(tmp)[3:5] = c("converged","CS","PIP")
tmp$method = 1
tmp2 = tmp
tmp = as.data.frame(rbind(tmp1,tmp2))

tmp$cs = apply(tmp[,c("sentinel","method","CS")],1,paste,collapse="-")

length(unique(tmp$sentinel))
length(unique(tmp$cs))
length(unique(subset(tmp,method=="1")$cs))
length(unique(subset(tmp,method=="10")$cs))

tab = table(tmp$cs)
tab = tab[tab==1]
length(tab)
tmp.cs1 = subset(tmp,cs %in% names(tab))
length(unique(tmp.cs1$snp))
length(unique(tmp.cs1$sentinel))
common.sub = subset(common,snp_id %in% tmp.cs1$snp)
nrow(common.sub)

cell = "Cluster24"
dataset = "corces_2020"
common.sub[,paste0("peak_overlap.",cell,".",dataset)]
idx = common.sub[,paste0("peak_overlap.",cell,".",dataset)]
common.sub[idx,paste0("abs_logfc.mean.",cell,".",dataset)]
common.sub[idx,paste0("abs_logfc.mean.pval.",cell,".",dataset)]








