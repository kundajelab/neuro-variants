# #############################################
# 
# # functions
# trim_cor_mat = function(cor.mat.tmp,thres=0.7) { 
#   # thres = 0.7
#   DONE=FALSE; 
#   while (!DONE) {
#     val = apply(cor.mat.tmp,2,function(x){sum(x > thres)})
#     for (i in seq(length(val),1)) {
#       if (val[i] != 1) {
#         cor.mat.tmp = cor.mat.tmp[-i,-i]
#         break
#       }
#       if (i==1) {DONE=TRUE}
#     }
#   }
#   return(cor.mat.tmp)
# }
# ########################################################
# 
# # trim to less-correlated cell types for cell-type-specificity analyses
# dataset_lst = c("corces_2020","trevino_2021")
# keep = c(); thres_val=0.6
# for (dataset in dataset_lst) {
#   cell_types = colnames(df)[grep("peak_overlap.",colnames(df))]
#   cell_types <- gsub("peak_overlap\\.", "", cell_types)
#   cell_types = grep(dataset,cell_types,value = TRUE)
#   unique_cell_types <- unique(cell_types)
#   
#   cols = paste0("abs_logfc.mean.",unique_cell_types)
#   cor.mat = cor(df[,cols],use = "na.or.complete")
#   colnames(cor.mat)=unique_cell_types
#   rownames(cor.mat)=unique_cell_types
#   cols = unique_cell_types
#   
#   cor.mat.tmp=cor.mat
#   cor.mat.tmp = trim_cor_mat(cor.mat.tmp,thres = thres_val)
#   keep = c(keep,rownames(cor.mat.tmp))
#   
# }
# 
# cell_types = colnames(df)[grep("peak_overlap.",colnames(df))]
# cell_types <- gsub("peak_overlap\\.", "", cell_types)
# cell_types = grep(paste(dataset_lst,collapse="|"),cell_types,value = TRUE)
# cols = paste0("abs_logfc.mean.",unique_cell_types)
# cor.mat = cor(df[,cols],use = "na.or.complete")
# colnames(cor.mat)=unique_cell_types
# rownames(cor.mat)=unique_cell_types
# cor.mat.sub = cor.mat[keep,keep]
# f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/cor_mat.txt"
# saveRDS(cor.mat,f.out)
# f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/cor_mat.trimmed.txt"
# saveRDS(cor.mat.sub,f.out)
# 
# # create new scores based on uncorrelated brain cell types
# keep.full = keep
# f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.ab_fb.txt"
# fwrite(as.data.frame(x=keep.full),f.out,row.names = F,col.names = F,quote = F)

f = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.ab_fb.txt"
keep.full = fread(f,data.table = F,stringsAsFactors = F,header = F)
keep.full = keep.full[,1]

# adult + fetal
keep = keep.full
df$cbp_min_pval_afb = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,min,na.rm=T)
cols <- paste0("abs_logfc.mean.", keep)
df$cbp_max_score_afb = apply(df[,cols],1,max,na.rm=T)
df$num_peaks_afb = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_afb = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})

# adult
keep = grep("corces_2020",keep.full,value = TRUE)
df$cbp_min_pval_ab = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,min,na.rm=T)
cols <- paste0("abs_logfc.mean.", keep)
df$cbp_max_score_ab = apply(df[,cols],1,max,na.rm=T)
df$num_peaks_ab = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_ab = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})

# fetal
keep = grep("trevino_2021",keep.full,value = TRUE)
df$cbp_min_pval_fb = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,min,na.rm=T)
cols <- paste0("abs_logfc.mean.", keep)
df$cbp_max_score_fb = apply(df[,cols],1,max,na.rm=T)
df$num_peaks_fb = apply(df[,paste0("peak_overlap.",keep)],1,sum,na.rm=T)
df$num_cbp_fb = apply(df[,paste0("abs_logfc.mean.pval.",keep)],1,function(x) {sum(x<0.01,na.rm = T)})

f.out=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot3.txt")
fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

