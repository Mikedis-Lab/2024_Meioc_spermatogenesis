#Maria Mikedis
#maria.mikedis@cchmc.org


#This script will 1) Import a seurat object (seu); 2) convert it into a monocle object (cds); and 3) save the monocle object


library(ggplot2)
library(data.table)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
theme_set(theme_bw())
library(dplyr)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(gplots)
library(RCurl)
library(monocle3)

###Set base directory###
base_dir <- "MM337_CellRanger_quant/filtered_feature_bc_matrix/integrated_results_WT.MeiocKO.only"
seu_object <- "integrated_data_processed_labeled_cellcycle.rds"

####Read in integrated Seurat object###
setwd(base_dir)
seu <- readRDS(file=seu_object)
DefaultAssay(seu) <- "RNA"
seu$celltype <- Idents(seu)

## subset to germ cell celltypes only
seu <- subset(seu, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut"
                                  )))

### convert Seurat object to a monocle object; https://github.com/cole-trapnell-lab/monocle3/issues/262

library("org.Mm.eg.db")
gene_symbol <- as.list(org.Mm.egSYMBOL)
corrected_data <- GetAssayData(seu, assay = "integrated", slot = "data")
raw_count_data <- GetAssayData(seu, assay = "RNA", slot = "counts")
class(raw_count_data)

cells_info <- seu@meta.data

gene_name <- gene_symbol[rownames(raw_count_data)]  
gene_name <- sapply(gene_name, function(x) x[[1]][1])

#preparing cds
gene_name <- ifelse(is.na(gene_name), names(gene_name), gene_name)
gene_short_name <- gene_name
gene_id <- rownames(raw_count_data)


genes_info <- cbind(gene_id, gene_short_name)
genes_info <- as.data.frame(genes_info)
rownames(genes_info) <- rownames(raw_count_data)

cds <- new_cell_data_set(expression_data = raw_count_data,
                         cell_metadata = cells_info,
                         gene_metadata = genes_info)

#double check that data looks ok
cds

#replace monocle clusters with seurat
cds@clusters$UMAP$partitions <- seu@meta.data$seurat_clusters
names(cds@clusters$UMAP$partitions) <- rownames(seu@meta.data)
cds@clusters$UMAP$clusters <- seu@meta.data$seurat_clusters
names(cds@clusters$UMAP$clusters) <- rownames(seu@meta.data)

#replace monocle dimensions with seurat
#cds@reducedDims$TSNE <-  seu@reductions$tsne@cell.embeddings
#cds@reducedDims$UMAP <-  seu@reductions$umap@cell.embeddings



### pseudotime
## Step 1: Normalize and pre-process the data ## https://cole-trapnell-lab.github.io/monocle3/docs/starting/
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells - browser window will pop up for selecting the starting branch of the trajectory
cds <- order_cells(cds, verbose=TRUE)



#Make new directory to hold violin plots
if (!dir.exists("pseudotime_monocle_GC.only")) {
  system(sprintf("mkdir pseudotime_monocle_GC.only"))
}


###SAVE LABELED RDS FILE###
saveRDS(cds, file = "pseudotime_monocle_GC.only/integrated_data_processed_labeled_monocle_WT.MeiocKO.only.rds")



