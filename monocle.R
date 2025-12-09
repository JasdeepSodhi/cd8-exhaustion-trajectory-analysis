library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(ggplot2)

#Input and Process Seurat to Monocle Object
data <- readRDS("CD8_DoubletRemoved_WholeAnno.RDS")
md <- data@meta.data
md <- md[, !duplicated(names(md)), drop = FALSE]      
names(md) <- make.unique(names(md))                 
data@meta.data <- md

cds <- as.cell_data_set(data, assay = "RNA", slot = "data")
cds <- estimate_size_factors(cds)
cd <- as.data.frame(colData(cds))
cd <- cd[, !duplicated(names(cd)), drop = FALSE]    
names(cd) <- make.unique(names(cd))                 
colData(cds) <- S4Vectors::DataFrame(cd)

fData(cds)$gene_short_name <- rownames(fData(cds))

recreate_partitions <- c(rep(1,length(cds@colData@rownames)))
names(recreate_partitions) <- cds@colData@rownames
recreate_partitions <- as.factor(recreate_partitions)

cds@clusters$UMAP$partitions <- recreate_partitions

list_clusters <- data@active.ident
cds@clusters$UMAP$clusters <- list_clusters

cds@int_colData@listData$reducedDims$UMAP <- data@reductions$umap@cell.embeddings

cds <- learn_graph(cds)
names(pData(cds)) <- make.unique(names(pData(cds)))


#Colors scale
colors <- c(
  "#FFEDFA",  # light pink
  "#FFB8E0",  # medium pink
  "#EC7FA9",  # pink-magenta
  "#BE5985",  # deeper magenta
  "#BFECFF",  # light sky blue
  "#CDC1FF",  # lavender
  "#FFF6E3",  # light cream
  "#FFE29A",  # light yellow
  "#FFCC80",  # peach
  "#FABEAF",  # salmon
  "#C8E7A4",  # soft green
  "#A2D6B4",  # mint green
  "#A6D8FF",  # baby blue
  "#C2B9F2",  # pastel indigo
  "#E4C1F9",  # light purple
  "#D3BCC0",  # dusty rose
  "#B5DED2",  # turquoise tint
  "#FFE5B4",  # melon
  "#EED5F5",  # lilac blush
  "#C1E1C1"   # soft sage
)

#FULL UMAP
a <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_cell_groups = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           show_trajectory_graph = FALSE) + 
  theme(legend.position = "right") +
  scale_color_manual(values = colors)


#Order cells (pseudotime) with CD8_7 cluster being initial
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == "CD8_7"]))
     
#Pseudotime full UMAP w.o trajectory
b <- plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_cell_groups = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           group_label_size = 10,
           show_trajectory_graph = FALSE)

#Pseudotime full UMAP w/ trajectory
c <- plot_cells(cds,
                     color_cells_by = 'pseudotime',
                     label_groups_by_cluster = FALSE,
                     label_cell_groups = FALSE,
                     label_branch_points = FALSE,
                     label_leaves = FALSE,
                     label_roots = FALSE,
                     group_label_size = 10,
                     show_trajectory_graph = TRUE,
                trajectory_graph_color = "black",     
                trajectory_graph_segment_size = 0.8) 
  


keep_clust <- c("CD8_7", "CD8_8", "CD8_4", "CD8_9", "CD8_3")
genes_of_interest <- c("CTLA4", "TCF7", "Znf683", "IFNG", "TOX", "GZMK", "LAG3", 
                       "PDCD1", "GZMH;GZMB")


#Keep only genes and clusters that are relevant  
cds_genes_filtered <- cds[rowData(cds)$gene_short_name %in% genes_of_interest, 
             colData(cds)$anno %in% keep_clust]
cds_genes_filtered <- order_cells(cds_genes_filtered, root_cells = colnames(cds_genes_filtered[,clusters(cds_genes_filtered) == "CD8_7"]))

#Genes in pseudotime
pseudotime_genes <- plot_genes_in_pseudotime(
  cds_genes_filtered,
  color_cells_by = "ident",
  min_expr = 0.3
 # trend_formula = "~ splines::ns(pseudotime, df = 5)"
) + 
  scale_color_manual(values = colors)

pseudotime_genestwo <- plot_genes_in_pseudotime(
  cds_genes_filtered,
  color_cells_by = "pseudotime",
  min_expr = 0.3
)

#Umap with genes of interest
genes_umap <- plot_cells(cds_genes_filtered,
           genes = genes_of_interest,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           min_expr = 0.1)

#Umaps for populations present in specific tissue sites
floating_umap1 <- plot_cells(
  cds,
  color_cells_by = "anno",
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  show_trajectory_graph = FALSE,   
  label_roots = FALSE
) +
  facet_wrap(~Tissue_Type, ncol = 2) +
  scale_color_manual(values = colors, breaks = levels(cds$anno)) +
  coord_equal() +
  labs(x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

floating_umap1$layers[[1]]$aes_params$alpha <- 0.3   

floating_umap1



