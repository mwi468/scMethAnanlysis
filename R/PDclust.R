#' @title  Clustering Function
#'
#' @description Clustering FNXN
#' @param  input containing c("Sample_ID","chr","start","end", "meth") columns
#' @import dplyr tidyverse devtools openxlsx PDclust parallel gridExtra pheatmap


Cluster_Methylation <- function(Methylation_table) {

dat <- Methylation_table[ , 1:5]
colnames(dat) <- c("Sample_ID","chr","start","end", "meth")
dat <- dat %>% filter(!grepl("GL|X|Y|MT", chr))%>%
  mutate(chr = paste("chr", chr, sep=""))


X<-split(dat, dat$Sample_ID)

X <- lapply(X, function(x) x[!(names(x) %in% "Sample_ID")])

PD <- create_pairwise_master(X)
PDmatrix<- convert_to_dissimilarity_matrix(PD)

summary(PD$num_cpgs)

## Generate clusters with PDclust.
PD_clusters <- cluster_dissimilarity(PDmatrix, num_clusters = 6)

## Determine cluster assignments.
PD_clust <- as.data.frame(PD_clusters$cluster_assignments)
PD_clust$sample_barcode <- rownames(PD_clust)


heatmap_pallete <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21)
pheatmap(PDmatrix,
         cluster_rows = PD_clusters$hclust_obj,
         cluster_cols = PD_clusters$hclust_obj,
         treeheight_row = 0,
         border_color = NA,
         color = heatmap_pallete,
         show_colnames = F,
         show_rownames = F,
         annotation_col = PD_clusters$cluster_assignments)

PDc_df <- visualize_clusters(PDmatrix, cluster_labels = PD_clusters$cluster_assignments)
PDc_df_revised <- as.data.frame(PDc_df)


plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background =element_rect(fill="white"))

ggplot(PDc_df_revised , aes(V1, V2, color = cluster)) +
  geom_point(alpha=0.8) +
  labs(color = "cluster", x="MDS Dimension 1", y="MDS Dimension 2") +
  plot_theme +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


}
