# Scree Plot
scree_df <- data.frame(PC = 1:length(pca_res$sdev), 
                       Variance = pca_res$sdev^2 / sum(pca_res$sdev^2) * 100)

ggplot(scree_df, aes(x = PC, y = Variance)) +
   geom_point(size = 5) +
   geom_line() +
   labs(title = "", x = "Principal Component", y = "Variance Explained (%)") +
   theme_minimal(base_size = 15) +
   theme(plot.title = element_text(hjust = 0.5))


ggplot(scree_df, aes(x = PC, y = Variance)) +
   geom_point(size = 5, color = "blue", alpha = 0.7) +  
   geom_line(color = "darkblue", size = 1) +  
   labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
   theme_minimal(base_size = 15) +
   theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  
         axis.title = element_text(size = 14, face = "bold"), 
         axis.text = element_text(size = 12),  
         panel.grid.major = element_line(color = "gray", linetype = "dashed"), 
         panel.grid.minor = element_blank(), 
         panel.border = element_rect(color = "black", fill = NA)) 

# Save as PDF
ggsave("scree_plot.pdf", plot = last_plot(), device = "pdf", width = 10, height = 8, dpi = 300)

# Save as JPEG
ggsave("scree_plot.jpg", plot = last_plot(), device = "jpeg", width = 10, height = 8, dpi = 300)

# Save as PNG
ggsave("scree_plot.png", plot = last_plot(), device = "png", width = 10, height = 8, dpi = 300)

# Save as TIFF
ggsave("scree_plot.tiff", plot = last_plot(), device = "tiff", width = 10, height = 8, dpi = 300)
