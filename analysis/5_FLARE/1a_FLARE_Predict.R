# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

# load
library(data.table)
library(glmnet)
# library(ranger)
# library(parallel)

################################################################################

# model list
model_lst = c("baseline",
              "baseline + fb peaks",
              "baseline + fb peaks + cbp",
              "baseline + brain peaks + cbp",
              "complete",
              "baseline + heart peaks + cbp",
              "baseline + ab peaks + cbp"
)

#######################

# input arguments
variantSet="asd"
bias="K562_bias"
i=3  

initial_data_load = function(variantSet) {
  # read dataframe
  cat("Reading data...\n")
  f=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot2.txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  
  # Filter SNPs with missing PhyloP values
  cat("Filter SNPs with missing PhyloP values... \n")
  ind = !is.na(df$phylop)
  df = df[ind,]
  
  cat("Performing misc. filtering... \n")
  
  # Use log10 distance to nearest TSS in the model:
  df$gene_distance_1.log10 = log10(df$gene_distance_1+1)
  
  # Create dummy variables for the categorical variable
  ind = !grepl("splice",df$Consequence)
  cat("Removing ",nrow(df) - sum(ind)," splice-related SNPs...\n")
  df = df[ind,]
  
  # return data
  return(df)
}

make_FLARE_input_data = function(variantSet,i,loaded_data=NULL) {
  
  if (is.null(loaded_data)) {
    df = initial_data_load(variantSet)
  } else {
    df = loaded_data
  }
  
  # select model
  model = model_lst[i]
  print(paste0("Model ",i,"/",length(model_lst),": ",model))
  print("...")
  
  baseline_cols = c("s_het_1","gene_distance_1.log10")
  
  if (model=="baseline") {
    cols = baseline_cols
  } else if (model=="baseline + fb peaks") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    summary_cols = c("num_peaks_fb")
    # cols = c(baseline_cols,peak_cols,summary_cols)
    cols = c(baseline_cols,peak_cols)
  } else if (model=="baseline + fb peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("trevino_2021|domcke_2020",cbp_cols,value = TRUE)
    cbp_cols <- grep("heart", cbp_cols, invert = TRUE, value = TRUE)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + ab peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("corces_2020",peak_cols,value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("corces_2020",cbp_cols,value = TRUE)
    summary_cols = c("num_cbp_ab","num_peaks_ab","cbp_max_score_ab")
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + brain peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("corces_2020|trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("corces_2020|trevino_2021|domcke_2020",cbp_cols,value = TRUE)
    cbp_cols <- grep("heart", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp_|num_peaks_|cbp_max_score_",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="complete") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp|num_peaks|cbp_max_score",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + heart peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("domcke_2020|ameen_2022|encode_2024",peak_cols,value = TRUE)
    peak_cols <- grep("brain", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("domcke_2020|ameen_2022|encode_2024",cbp_cols,value = TRUE)
    cbp_cols <- grep("brain", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp_|num_peaks_|cbp_max_score_",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  }
  
  ##############################################################################
  # Interaction terms:
  cat("Creating interaction terms...\n")
  
  # This line excludes all columns in df whose names start with int_. 
  # The code resets df in each iteration to a version of df without those columns.
  df = df[,!(colnames(df) %in% grep("^int_",colnames(df),value=TRUE))] # reset in a for loop
  if (!(model %in% c("baseline","baseline + fb peaks"))) {
    for (k in 1:length(peak_cols)) {
      # cat(k,"/",length(peak_cols),"\n",sep = '')
      df[,paste0("int_",k)] = df[,cbp_cols[k]] * df[,peak_cols[k]]
    }
    cols = c(cols,paste0("int_",1:length(peak_cols)))
  }
  
  return(list(df,cols))
}

##############################################################################

FLARE_Predict = function(variantSet,i,loaded_data=NULL) {
  
  # initialize:
  model = model_lst[i]
  out = make_FLARE_input_data(variantSet,i,loaded_data)
  df = out[[1]]
  cols = out[[2]]
  rm(out)
  predictions_df.all = list(); cor_result = list()
  
  for (chrNum in 1:22) {
    cat("Making FLARE predictions using LOCO schema (chr ",chrNum,")...\n",sep = '')
    # Subset to chr of interest
    ind_chr_include = df$chr==paste0("chr",chrNum)
    x = as.matrix(df[ind_chr_include,cols])
    # Load model
    f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/models/",i,".",model,".chr",chrNum,".lasso.rds")
    final_mod = readRDS(f)
    # Making FLARE predictions
    predictions_lasso <- as.numeric(predict(final_mod, newx = x)[,1])
    # Store predictions
    predictions_df = data.frame(snp_id = df[ind_chr_include,"snp_id"],predictions_lasso)
    # Store results
    predictions_df.all[[chrNum]] = predictions_df
  }
  # Save predictions
  predictions_df = as.data.frame(do.call(rbind,predictions_df.all))
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/models/",i,".",model,".predictions_",variantSet,".lasso.txt")
  fwrite(predictions_df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  return(predictions_df$predictions_lasso[match(df$snp_id, predictions_df$snp_id)])
  
}

################################################################################
################################################################################
################################################################################

variantSet = "common"
for (variantSet in c("asd","common","rare")) {
  modelNum = 1
  df = initial_data_load(variantSet)
  df$FLARE_baseline <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 2
  df$FLARE_fb_peaks <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 3
  df$FLARE_fb <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 4
  df$FLARE_brain <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 5
  df$FLARE_all <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 6
  df$FLARE_heart <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  modelNum = 7
  df$FLARE_ab <- FLARE_Predict(variantSet,modelNum,loaded_data=df)
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/",variantSet,".FLARE.txt")
  fwrite(df[,c("snp_id","phylop","FLARE_baseline","FLARE_fb_peaks","FLARE_fb","FLARE_brain","FLARE_all","FLARE_heart","FLARE_ab")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

cor(df[,c("phylop","FLARE_baseline","FLARE_fb_peaks","FLARE_fb","FLARE_heart","FLARE_ab")])

# Metrics for FLARE predictions:

print(cor(df$FLARE_fb_peaks,df$phylop))
print(cor(df$FLARE_fb,df$phylop))
tmp = make_FLARE_input_data("asd",3,loaded_data=df)[[1]]
tmp$accessible = apply(tmp[,grep("^int_",colnames(tmp))],1,sum) > 0
tmp = subset(tmp,accessible)
print(cor(tmp$FLARE_fb_peaks,tmp$phylop))
print(cor(tmp$FLARE_fb,tmp$phylop))
# chromBPnet versus peaks: 
# ~30% improvement genome-wide
# ~50% improvement at accessible variants



