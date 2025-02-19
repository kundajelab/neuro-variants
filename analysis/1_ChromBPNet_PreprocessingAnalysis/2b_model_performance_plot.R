# Local script:

# This script processes and visualizes chromBPnet 
# model performance across different cellular contexts, 
# generating a publication-ready bar plot.
  
################################################################################
  
# Load necessary libraries
library(data.table) # For efficient data manipulation
library(ggplot2)    # For data visualization

# Read the main dataset containing chromBPnet model performance metrics
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/data/cbp/model_performance.tsv",
           data.table = F, stringsAsFactors = F)

# Alternative directory path (uncomment if using an alternate setup)
# df = fread("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/cbp/model_performance.tsv",
#            data.table = F, stringsAsFactors = F)

# Rename the first column for clarity
colnames(df)[1] = "cell"

# Create a new column combining development stage and organ type
df$condition = paste(df$dev, df$organ)

# Read and merge cluster descriptions from external label file
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/neuro-variants/output/data/etc/adult_cluster_names.csv"
labels = fread(f, data.table = F, stringsAsFactors = F)
df = merge(df, labels[, c("Cluster", "Cluster_Description")], by.x = "cell", by.y = "Cluster", all.x = TRUE)

# Initialize and populate cell type annotations
df$cell_v2 = NA
df$cell_v2[!is.na(df$Cluster_Description)] = df$Cluster_Description[!is.na(df$Cluster_Description)]

# Merge additional cell annotations from another external label file
f = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/neuro-variants/output/data/etc/cell_cluster_annot.csv"
labels = fread(f, data.table = F, stringsAsFactors = F)
df = merge(df, labels[, c("Cluster ID", "Name (long)")], by.x = "cell", by.y = "Cluster ID", all.x = TRUE)

# Update `cell_v2` annotations if new labels are available
df$cell_v2[!is.na(df$`Name (long)`)] = df$`Name (long)`[!is.na(df$`Name (long)`)]

# Remove redundant columns
df = df[, !(colnames(df) %in% c("Name (long)", "Cluster_Description"))]

# Set `cell_v2` to `cell` for rows with no annotations
df$cell_v2[is.na(df$cell_v2)] = df$cell[is.na(df$cell_v2)]

# Identify neuron-related cell types
df$neuron = grepl("neuron", df$cell_v2, ignore.case = TRUE)

# Assign human-readable context labels for different datasets and conditions
df$context = NA
df$context[df$dataset == "domcke_2020" & df$neuron & df$condition == "fetal brain"] = "Fetal brain neurons (Domcke)"
df$context[df$dataset == "domcke_2020" & !df$neuron & df$condition == "fetal brain"] = "Fetal brain non-neurons (Domcke)"
df$context[df$dataset == "domcke_2020" & df$condition == "fetal heart"] = "Fetal heart (Domcke)"
df$context[df$dataset == "corces_2020"] = "Adult brain (Corces)"
df$context[df$dataset == "encode_2024"] = "Adult heart (ENCODE)"
df$context[df$dataset == "ameen_2022"] = "Fetal heart (Ameen)"
df$context[df$dataset == "trevino_2021" & df$neuron & df$condition == "fetal brain"] = "Fetal brain neurons (Trevino)"
df$context[df$dataset == "trevino_2021" & !df$neuron & df$condition == "fetal brain"] = "Fetal brain non-neurons (Trevino)"

# Simplify context labels and set an ordered factor
df$context_v2 = sub("\\s*\\(.*", "", df$context)
df$context_v2 = factor(df$context_v2, levels = c("Fetal brain neurons", "Fetal brain non-neurons", "Fetal heart", "Adult heart", "Adult brain"))

# Identify fetal brain neurons as a separate category
df$fetal_neuron = df$context_v2 == "Fetal brain neurons"

# Save the full dataset for reference
df.full = df

# Aggregate Spearman correlation by cell and context for plotting
df = aggregate(peaks_spearmanr ~ cell + context_v2, df, mean)

# Define output file path for the figure
f.out = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/neuro-variants/output/figs/cbp_model_performance.pdf"

# Alternative output path (uncomment if using a different setup)
# f.out = "/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/cbp/analysis/plots/cbp_model_performance.pdf"

# Create and save the bar plot
pdf(f.out, width = 12, height = 3.5)
g = ggplot(df, aes(x = reorder(cell, -peaks_spearmanr, mean),
                   y = peaks_spearmanr, 
                   fill = context_v2)) +
  geom_bar(stat = 'identity') + # Bar plot showing average Spearman correlations
  scale_fill_manual(values = c("Fetal brain neurons" = "#E0CA70",
                               "Adult brain" = "#483FA3",
                               "Fetal brain non-neurons" = "#4D3B3B",
                               "Adult heart" = "#A34D3F",
                               "Fetal heart" = "#852222")) + # Custom colors
  ggpubr::theme_pubr() + # Publication-ready theme
  labs(x = "Performance of chromBPnet models across cellular contexts",
       y = "Spearman correlation (peaks)",
       fill = "") + # Axis labels and legend title
  theme(axis.text.x = element_blank()) + # Optional: Hide x-axis text
  # coord_flip() + # Optional: Rotate the plot
  # theme(axis.text.y = element_blank()) + # Optional: Hide y-axis text
  geom_point(data = df.full, aes(x = cell, y = peaks_spearmanr),
             position = position_jitter(width = 0.05), color = "black", size = rel(0.5),
             show.legend = FALSE) # Add jittered points for individual data

print(g)
dev.off()
