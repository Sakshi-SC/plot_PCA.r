# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)  # For better labeling
library(RColorBrewer)  # For more color options

setwd("/home/rohan/Documents/Sakshi/plots/pca_plot") #set working direcctory to save plots 

 dds <- DESeq(dds) #deg result

# Extract results for Test vs. Control
res <- results(dds, contrast = c("condition", "Test", "Control"), alpha = 0.05)

# Identify DEGs based on adjusted p-value and fold-change thresholds
degs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# Extract gene names for DEGs
deg_gene_names <- rownames(degs)

# Apply Variance-Stabilizing Transformation (VST) to the DESeq2 object
vsd <- vst(dds, blind = FALSE)

# Extract the VST-transformed expression matrix
expr_vst <- assay(vsd)

# Filter expression matrix to include only DEGs
expr_vst_degs <- expr_vst[rownames(expr_vst) %in% deg_gene_names, ]

# Remove zero-variance columns and rows 
zero_variance_cols <- apply(expr_vst_degs, 2, var) == 0
expr_vst_degs <- expr_vst_degs[, !zero_variance_cols]
zero_variance_rows <- apply(expr_vst_degs, 1, var) == 0
expr_vst_degs <- expr_vst_degs[!zero_variance_rows, ]

# Perform PCA on the filtered DEG data (transposed)
pca_res <- prcomp(t(expr_vst_degs), scale. = TRUE)

# Summary of PCA
print(summary(pca_res))

# data frame for plotting
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Sample = colnames(expr_vst_degs))

# Applying Variance-Stabilizing Transformation (VST) to the DESeq2 object
vsd <- vst(dds, blind = FALSE)

# Extract the VST-transformed expression matrix
expr_vst <- assay(vsd)

# Filter expression matrix to include only DEGs
expr_vst_degs <- expr_vst[rownames(expr_vst) %in% deg_gene_names, ]

# Remove zero-variance columns and rows 
zero_variance_cols <- apply(expr_vst_degs, 2, var) == 0
expr_vst_degs <- expr_vst_degs[, !zero_variance_cols]
zero_variance_rows <- apply(expr_vst_degs, 1, var) == 0
expr_vst_degs <- expr_vst_degs[!zero_variance_rows, ]

# Perform PCA on the filtered DEG data (transposed)
pca_res <- prcomp(t(expr_vst_degs), scale. = TRUE)

# Summary of PCA
print(summary(pca_res))

# data frame for plotting
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Sample = colnames(expr_vst_degs))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample, shape = Sample)) +
   geom_point(size = 19, alpha = 0.8) +  # this is size of point
   labs(title = "", 
        x = paste0("PC1: ", round(summary(pca_res)$importance[2,1] * 100, 1), "% variance"),
        y = paste0("PC2: ", round(summary(pca_res)$importance[2,2] * 100, 1), "% variance")) +
   theme_minimal(base_size = 20) +  # font size
   theme(legend.position = "right",  
         plot.title = element_text(hjust = 0.5), 
         axis.line = element_line(color = "black"), 
         axis.ticks = element_line(color = "black")) +  
   scale_color_brewer(palette = "Dark2") +  
   scale_shape_manual(values = rep(c(19, 19), each = 3)) +  
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) +  # Vertical line at PC1 = 0
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 1)    # Horizontal line at PC2 = 0

# Save as PDF
ggsave("pca_plot.pdf", plot = last_plot(), device = "pdf", width = 10, height = 8, dpi = 300)

# Save as JPEG
ggsave("pca_plot.jpg", plot = last_plot(), device = "jpeg", width = 10, height = 8, dpi = 300)

# Save as PNG
ggsave("pca_plot.png", plot = last_plot(), device = "png", width = 10, height = 8, dpi = 300)

# Save as TIFF
ggsave("pca_plot.tiff", plot = last_plot(), device = "tiff", width = 10, height = 8, dpi = 300)

