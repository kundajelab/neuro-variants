# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

# load
library(data.table)
library(glmnet)

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

FLARE_Training = function(variantSet,i,loaded_data=NULL) {
  
  # initialize:
  model = model_lst[i]
  out = make_FLARE_input_data(variantSet,i,loaded_data)
  df = out[[1]]
  cols = out[[2]]
  rm(out)
  
  for (chrNum in 1:22) {
    
    cat("Training FLARE using LOCO schema (chr ",chrNum,")...\n",sep = '')
    
    # Curate data:
    ind_chr_exclude = df$chr!=paste0("chr",chrNum)
    x = as.matrix(df[ind_chr_exclude,cols])
    y = df$phylop[ind_chr_exclude]
    
    ############################################################################
    # # Lasso
    set.seed(123)  # Set the seed for reproducibility:
    
    # Set peak, chrombpnet, and gene constraint info to have non-negative constraints.
    constrained_cols = unique(c(
      grep("abs_logfc.mean",colnames(x),value = TRUE),
      grep("cbp",colnames(x),value = TRUE),
      grep("peak",colnames(x),value = TRUE),
      grep("int_",colnames(x),value=TRUE),
      grep("s_het_1",colnames(x),value = TRUE)
    ))
    constrained_limits = rep(-Inf,ncol(x))
    constrained_limits[colnames(x) %in% constrained_cols] = 0
    
    # L1 loss with non-negative constraints - identify lambda values
    alphaUse = 1
    cv_mod <- cv.glmnet(x, y, alpha = alphaUse, lower.limits = constrained_limits, nfolds = 4)  # 4-fold cross-validation
    
    # Fit model with optimal lambda
    best_lambda <- cv_mod$lambda.min
    final_mod <- glmnet(x, y, alpha = alphaUse, lower.limits = constrained_limits, lambda = best_lambda)
    
    # Save model
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/models/",i,".",model,".chr",chrNum,".lasso.rds")
    saveRDS(final_mod,file=f.out)
  }
}

##############################################################################

df = initial_data_load("rare")
for (modelNum in 1:7) {
  print(modelNum)
  FLARE_Training("rare",modelNum,loaded_data=df)
}








