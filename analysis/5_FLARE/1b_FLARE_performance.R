library(data.table)
df1 = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/asd.FLARE.txt",data.table = F,stringsAsFactors = F)
df2 = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/rare.FLARE.txt",data.table = F,stringsAsFactors = F)
df3 = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/common.FLARE.txt",data.table = F,stringsAsFactors = F)

cor_result = cor.test(df1$phylop,df1[,"FLARE_fb"])
cor_result
df1.sub = subset(df1,FLARE_fb > -0.1)
cor_result = cor.test(df1.sub$phylop,df1.sub[,"FLARE_baseline"])
cor_result
cor_result = cor.test(df1.sub$phylop,df1.sub[,"FLARE_fb"])
cor_result
cor_result = cor.test(df1.sub$phylop,df1.sub[,"FLARE_fb_peaks"])
cor_result

subset(df1,)


predictor_lst = c("FLARE_baseline","FLARE_fb_peaks","FLARE_fb","FLARE_brain","FLARE_all","FLARE_heart","FLARE_ab")
res.lst = list(); k=0
for (i in 1:length(predictor_lst)) {
  predictor = predictor_lst[i]
  k = k+1
  cor_result = cor.test(df1$phylop,df1[,predictor])
  res.lst[[k]] = data.frame(variantSet="ASD",model=predictor,r=cor_result$estimate,l=cor_result$conf.int[1],h=cor_result$conf.int[2])
}
for (i in 1:length(predictor_lst)) {
  predictor = predictor_lst[i]
  k = k+1
  cor_result = cor.test(df2$phylop,df2[,predictor])
  res.lst[[k]] = data.frame(variantSet="rare",model=predictor,r=cor_result$estimate,l=cor_result$conf.int[1],h=cor_result$conf.int[2])
}
for (i in 1:length(predictor_lst)) {
  predictor = predictor_lst[i]
  k = k+1
  cor_result = cor.test(df3$phylop,df3[,predictor])
  res.lst[[k]] = data.frame(variantSet="common",model=predictor,r=cor_result$estimate,l=cor_result$conf.int[1],h=cor_result$conf.int[2])
}
res = as.data.frame(do.call(rbind,res.lst))

fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/flare_performance.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

res = list()

y = apply(df1[,-1],2,mean)
y2 = apply(df1[,-1],2,sd)
n = nrow(df1)
se = y2/sqrt(n)
res[[1]] = data.frame(variantSet="asd",metric=names(y),mu=y,se,row.names = NULL)
y = apply(df2[,-1],2,mean)
y2 = apply(df2[,-1],2,sd)
n = nrow(df2)
se = y2/sqrt(n)
res[[2]] = data.frame(variantSet="rare",metric=names(y),mu=y,se,row.names = NULL)
y = apply(df3[,-1],2,mean)
y2 = apply(df3[,-1],2,sd)
n = nrow(df3)
se = y2/sqrt(n)
res[[3]] = data.frame(variantSet="common",metric=names(y),mu=y,se,row.names = NULL)
res = as.data.frame(do.call(rbind,res))

fwrite(res,"/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/flare_mean.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



