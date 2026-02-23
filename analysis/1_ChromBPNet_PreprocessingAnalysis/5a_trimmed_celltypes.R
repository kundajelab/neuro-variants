# functions
trim_cor_mat = function(cor.mat.tmp,thres=0.7) {
  # thres = 0.7
  DONE=FALSE;
  while (!DONE) {
    val = apply(cor.mat.tmp,2,function(x){sum(x > thres)})
    for (i in seq(length(val),1)) {
      if (val[i] != 1) {
        cor.mat.tmp = cor.mat.tmp[-i,-i]
        break
      }
      if (i==1) {DONE=TRUE}
    }
  }
  return(cor.mat.tmp)
}
########################################################

# trim to less-correlated cell types for cell-type-specificity analyses
keep = c(); thres_val=0.7
# for (dataset in dataset_lst) {
for (i in 1:4) {
  if (i==1) {
    context="fetal_brain"
    unique_cell_types = fetal_brain
  } else if (i==2) {
    context="fetal_heart"
    unique_cell_types = fetal_heart
  } else if (i==3) {
    context="adult_brain"
    unique_cell_types = adult_brain
  } else if (i==4) {
    context="adult_heart"
    unique_cell_types = adult_heart
  }
  cols = paste0("abs_logfc.mean.",unique_cell_types)
  cor.mat = cor(rare[,cols],use = "pairwise.complete.obs")
  colnames(cor.mat)=unique_cell_types
  rownames(cor.mat)=unique_cell_types
  cols = unique_cell_types
  cor.mat.tmp=cor.mat
  cor.mat.tmp = trim_cor_mat(cor.mat.tmp,thres = thres_val)
  keep = rownames(cor.mat.tmp)
  print(length(keep))
  print(keep)
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/analysis/keep_cols.",context,".txt")
  fwrite(as.data.frame(x=keep),f.out,row.names = F,col.names = F,quote = F)
}
