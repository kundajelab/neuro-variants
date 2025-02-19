library(data.table)
traitName = "Alzheimers_Bellenguez_2022"
gwas = fread("/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/Alzheimers_Bellenguez_2022/Alzheimers_Bellenguez_2022.txt.gz",data.table = F,stringsAsFactors = F)
fm = fread(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/all/",traitName,".susie_all.txt"),data.table = F,stringsAsFactors = F)

fm.sub=subset(fm,sentinel=="chr5:180201150:G:A" & CS.L_1 > -1);fm.sub
fm.sub=subset(fm,sentinel=="chr5:180201150:G:A")
gwas.sub = subset(gwas,variant_alternate_id %in% fm.sub$snp)

f.out = "/oak/stanford/groups/smontgom/amarder/tmp/alzheimers_locus.gwas.txt"
fwrite(gwas.sub,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
f.out = "/oak/stanford/groups/smontgom/amarder/tmp/alzheimers_locus.finemap.txt"
fwrite(fm.sub,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

################################################################

gwas.sub = fread("~/Downloads/alzheimers_locus.gwas.txt",data.table = F,stringsAsFactors = F)
fm.sub = fread("~/Downloads/alzheimers_locus.finemap.txt",data.table = F,stringsAsFactors = F)
head(fm.sub)
head(gwas.sub)

fm.sub$pipmax = apply(fm.sub[,c("PIP.L_10","PIP.L_1")],1,max)
df.mg = merge(gwas.sub[,c("variant_alternate_id","base_pair_location","p_value")],fm.sub[,c("snp","pipmax","PIP.L_10","PIP.L_1","CS.L_1","CS.L_10")],by.x="variant_alternate_id",by.y="snp")
df.mg[(df.mg$pipmax > 0.2),]
df.mg.sub = subset(df.mg,base_pair_location > (180201150-50000) &
                     base_pair_location < (180201150+50000))
df.mg.sub[(df.mg.sub$pipmax > 0.2),]

g1 = ggplot(df.mg.sub,aes(x=base_pair_location,y=-log10(p_value),color=as.factor(CS.L_1))) +
  geom_point(alpha=0.4) +
  geom_hline(yintercept = -log10(5e-8),col='red',lty='dashed') +
  labs(x="Position",y= bquote(-log[10]~italic(P))) +
  theme_bw() +
  scale_color_manual(values=c("black","orange")) +
  guides(col="none") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8));g1
g2 = ggplot(df.mg.sub,aes(x=base_pair_location,y=pipmax,color=as.factor(CS.L_1))) +
  geom_point(alpha=0.4) +
  labs(x="Position",y= "PIP") +
  theme_bw() +
  scale_color_manual(values=c("black","orange")) +
  guides(col="none") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8));g2

pdf("~/Downloads/alzheimers_locus.gwas.pdf",width = 5*2,height=0.8*2)
print(g1)
dev.off()
pdf("~/Downloads/alzheimers_locus.finemap.pdf",width = 5*2,height=0.6*2)
print(g2)
dev.off()

dprint(g2)

