# Cluster script
# conda activate r

# Load the data.table library for efficient data manipulation
library(data.table)
library(dplyr)

# Initialize an empty list to hold data frames from each dataset
df.lst = list()

# Loop through datasets to load model performance data
for (dataset in c("ameen_2022", "corces_2020", "domcke_2020", "encode_2024", "trevino_2021")) {
  # Construct the file path for the current dataset
  f = paste0("/oak/stanford/groups/akundaje/projects/neuro-variants/model_performance/", dataset, "/model_performance.tsv")
  
  # Load the dataset into a data frame
  df = fread(f, data.table = F, stringsAsFactors = F)
  
  # Add a column to identify the dataset
  df$dataset = dataset
  
  # Store the data frame in the list
  df.lst[[dataset]] = df
}

# Combine all dataset-specific data frames into one large data frame
df = as.data.frame(do.call(rbind, df.lst))

# Initialize a new column to classify organ type
df$organ = "heart"

# Classify organ as "brain" for specific datasets or cell types
i = (df$dataset %in% c("corces_2020", "trevino_2021")) | (df$dataset == "domcke_2020" & grepl("brain", df$celltype))
df$organ[i] = "brain"

# Initialize a new column to classify developmental stage
df$dev = "fetal"

# Classify developmental stage as "adult" for specific datasets
i = (df$dataset %in% c("corces_2020", "encode_2024"))
df$dev[i] = "adult"

# Calculate mean performance metrics (Pearson r, Spearman r, median JSD) grouped by cell type and dataset
df <- df %>%
  group_by(celltype, dataset) %>%
  mutate(
    mean_peaks_pearsonr = mean(peaks_pearsonr, na.rm = TRUE),
    mean_peaks_spearmanr = mean(peaks_spearmanr, na.rm = TRUE),
    mean_peaks_median_jsd = mean(peaks_median_jsd, na.rm = TRUE)
  ) %>%
  ungroup() %>% as.data.frame()


# remove outlier that are > 4 SD outside typical model performance using spearman r
# df$outlier_score = abs(scale(df$mean_peaks_spearmanr))
# df[order(df$outlier_score,decreasing = T),][1:3,]
tmp = df[!duplicated(paste(df$celltype,df$dataset)),]
tmp = tmp[abs(scale(tmp$mean_peaks_spearmanr)) > 4,]
outlier_remove = df[(df$celltype %in% tmp$celltype),]
# df = df[abs(scale(df$peaks_spearmanr)) < 5,]
df = df[!(df$celltype %in% outlier_remove$celltype),]
# "LV.adipocyte.In" and "LV.mast_cell.In" removed

# save
f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/model_performance.tsv"
fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/cbp/model_performance.outliers.tsv"
fwrite(outlier_remove,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# print some summary stats:
tmp = aggregate(peaks_pearsonr ~ celltype,df,mean)
median(tmp$peaks_pearsonr)
mean(tmp$peaks_pearsonr)
tmp = aggregate(peaks_spearmanr ~ celltype,df,mean) # spearman might be the most "fair" and interpretable benchmark
mean(tmp$peaks_spearmanr)
median(tmp$peaks_spearmanr)
range(tmp$peaks_spearmanr)
tmp = aggregate(peaks_median_jsd ~ celltype,df,mean)
mean(tmp$peaks_median_jsd)



