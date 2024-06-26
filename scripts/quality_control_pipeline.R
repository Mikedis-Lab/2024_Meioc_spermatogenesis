#Maria Mikedis
#maria.mikedis@cchmc.org


#GENERATES:
#1) Basic QC violinplots (% ribo, % mito, features, counts, % doublets)
#2) Table showing # cells for each sample (with different degrees of filtering)
#3) Plot showing distribution of doublets


###SET PARAMETERS###
base_directory <- "MM337_CellRanger_quant/filtered_feature_bc_matrix/"   #directory with files matrix.mx.gz, features.tsv.gz, and barcodes.tsv.gz, one per sample
sample_names <- c("L16_2807", "L16_2808", "L16_3018", "L16_3020")      #Names of the relevant sub-folders in bus_output 
metadata <- "genome_files/Mus_musculus_Gencode.vM23.basic/gencode.vM23.basic_gcna_txo_annotable.txt"  #tab-separated gene-level metadata file containing, at minimum, column with gene_name and column with gene_type (i.e."protein-coding")

###IMPORT PACKAGES###
library(magrittr)
library(DropletUtils)
library(reticulate)
library(BUSpaRse)
library(ggplot2)
library(data.table)
library(Seurat)
library(Matrix)
theme_set(theme_bw())
library(dplyr)
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(reshape)
library(tidyverse)
library(scales)
library(RCurl)


memory.limit(27000) ### to provide enough memory to run on cluster; can run in R through tak without it


###READ CELLRANGER OUTPUT INTO SEURAT OBJECT
#https://github.com/hbctraining/scRNA-seq/blob/master/lessons/03_SC_quality_control-setup.md


# Create Seurat object with all samples
#for (file in sample_names){
#        seurat_data <- Read10X(data.dir = paste(base_directory, file, sep=""))
#        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
#                                         min.features = 100, 
#                                         project = file)
#        assign(file, seurat_obj)
#}

#Check each sample

#head(L16_2807@meta.data)
#head(L16_2808@meta.data)


#as.data.frame(seurat_obj@assays$RNA@counts[1:10, 1:2])



####READ IN METADATA###
meta <- read.csv(metadata, sep="\t")
meta <- distinct(meta, gene_name, .keep_all = TRUE)
rownames(meta) <- meta$gene_name

genes_of_interest <- c()

for (i in 1:length(rownames(meta))) {
  gene <- rownames(meta)[i]
  if (meta[gene,"gene_type"]=="protein_coding") {
    genes_of_interest <- c(genes_of_interest,gene)
  }
  else if (gene %in% c("Xist")) {
    genes_of_interest <- c(genes_of_interest,gene)
  }
}


###PROCESS EACH SAMPLE INDIVIDUALLY###
setwd(base_directory)

#Get list of processed Seurat objects for integration
object_list <- c()

for (sample in sample_names) {
  
  #Read in file
  #directory <- sprintf("bus_output/%s/genecount", sample)
  res_mat <- Read10X(data.dir = paste(base_directory, sample, sep=""))
  
  #Filter to protein-coding genes only
  res_mat <- res_mat[row.names(res_mat) %in% genes_of_interest,]
  
  #Make object
  object <- CreateSeuratObject(counts=res_mat, min.cells = 3, min.features=200, project=sample) 
  
  #Make violin plot and counts plot for unfiltered data
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^Rp[ls]")
  object[["percent.germcell"]] <- PercentageFeatureSet(object, features = c("Dazl","Ddx4"))

  
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object<- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  
  #Add column with doublet scores
  sce <- as.SingleCellExperiment(object)
  sce <- cxds(sce,retRes = TRUE)  #coexpression
  sce <- bcds(sce,retRes = TRUE,verb=TRUE) #artificial doublets
  sce <- cxds_bcds_hybrid(sce)
  
  object@meta.data$hybrid_score <- sce$hybrid_score
  object@meta.data$bcds_score <- sce$bcds_score
  object@meta.data$cxds_score <- sce$cxds_score
  
  object_list <- c(object_list, object)
   
}

###INTEGRATE SAMPLES###
seu.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:25, reference = c(1)) #takes around 15 min for 4 samples
seu <- IntegrateData(anchorset = seu.anchors, dims = 1:25)

DefaultAssay(seu) <- "integrated"

###SAVE FILE###
#saveRDS(seu, file = "quality_control/integrated_data_unprocessed.rds")

###QC VIOLIN PLOTS###
if (!dir.exists("quality_control")) {
  system(sprintf("mkdir quality_control"))
}

#
png("quality_control/features_versus_mito.png", width = 1480, height = 740, units = "px")
FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=0.5) + 
  geom_hline(yintercept=15, linetype="dashed", color = "blue", size=0.5) + 
  geom_vline(xintercept=1000, linetype="dashed", color = "black", size=0.5) + 
  geom_vline(xintercept=500, linetype="dashed", color = "blue", size=0.5) 
dev.off()

#Trim the top 1%:  Otherwise outliers shift scale
top_1_features <- quantile(seu@meta.data$nFeature_RNA, probs = c(0.99))
seu_features <- subset(seu, subset = nFeature_RNA < top_1_features)
png("quality_control/violin_integrated_features.png", width = 1480, height = 740, units = "px")
VlnPlot(seu_features, features = c("nFeature_RNA"),pt.size=0.01) + theme(legend.position = 'none')
dev.off()

#Trim the top 1%:  Otherwise outliers shift scale
top_1_counts <- quantile(seu@meta.data$nCount_RNA, probs = c(0.99))
seu_counts <- subset(seu, subset = nCount_RNA < top_1_counts)
png("quality_control/violin_integrated_counts.png", width = 1480, height = 740, units = "px")
VlnPlot(seu_counts, features = c("nCount_RNA"),pt.size=0.01) + theme(legend.position = 'none')
dev.off()

#Trim the top 5%:  Otherwise outliers shift scale
top_5_mito <- quantile(seu@meta.data$percent.mt, probs = c(0.95))
seu_mito <- subset(seu, subset = percent.mt < top_5_mito)
png("quality_control/violin_integrated_mito.png", width = 1480, height = 740, units = "px")
VlnPlot(seu_mito, features = c("percent.mt"),pt.size=0.01) + theme(legend.position = 'none')
dev.off()

png("quality_control/violin_integrated_ribo.png", width = 1480, height = 740, units = "px")
VlnPlot(seu, features = c("percent.ribo"),pt.size=0.01) + theme(legend.position = 'none')
dev.off()

png("quality_control/violin_integrated_doublets.png", width = 1480, height = 740, units = "px")
VlnPlot(seu, features = c("bcds_score"),pt.size=0.01) + theme(legend.position = 'none')
dev.off()


###CELL COUNT TABLE###

min_200 <- table(seu@meta.data$orig.ident)

seu_500 <- subset(seu, subset = nFeature_RNA > 500)
min_500 <- table(seu_500@meta.data$orig.ident)

seu_500_mito_10 <- subset(seu_500, subset = percent.mt <10)
min_500_mito_10 <- table(seu_500_mito_10@meta.data$orig.ident)

seu_1000 <- subset(seu, subset = nFeature_RNA > 1000)
min_1000 <- table(seu_1000@meta.data$orig.ident)

seu_1000_mito_10 <- subset(seu_1000, subset = percent.mt <10)
min_1000_mito_10  <- table(seu_1000_mito_10@meta.data$orig.ident)

total <- rbind(min_200, min_500, min_500_mito_10, min_1000, min_1000_mito_10)
write.table(total, file = "quality_control/cells_per_sample.csv", sep = ",", col.names=NA)

total <- melt(total)
total$X1 <- factor(total$X1, levels= c("min_200","min_500","min_500_mito_10","min_1000","min_1000_mito_10"))

png("quality_control/cells_per_sample.png", width = 1480, height = 740, units = "px")
ggplot(total, aes(factor(X2), value, fill = X1)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width=0.5) + 
  scale_fill_brewer(palette = "Set1") + xlab("") + ylab("Number of cells") + guides(fill=guide_legend(title="Threshold")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 20))
dev.off()

###Plot doublet scores as histograms and TSNEs###

#Plot TSNEs as grid
DefaultAssay(seu) <- "RNA"

make_doublet_tsne <- function(sample_name) {
  
  subset <- subset(seu, orig.ident==sample_name)
  
  sce <- as.SingleCellExperiment(subset)
  
  logcounts(sce) <- log1p(counts(sce))
  vrs <-apply(logcounts(sce),1,var)
  pc <-rpca(t(logcounts(sce)[order(vrs,decreasing=TRUE)[1:100],]))
  ts <-Rtsne(pc$x[,1:10],verb=FALSE,check_duplicates = FALSE)
  
  reducedDim(sce,"tsne") = ts$Y; rm(ts,vrs,pc)
  
  plotReducedDim(sce,"tsne",col="bcds_score") + ggtitle(sample_name)
  
}

list <-lapply(sample_names, make_doublet_tsne)

png("quality_control/doublet_tsne.png", width = 1500, height = 1000, units="px")
cowplot::plot_grid(plotlist = list)
dev.off()

#Plot histograms as grid
make_doublet_hist <- function(sample_name) {
  
  subset <- subset(seu, orig.ident==sample_name)
  
  qplot(subset@meta.data$bcds_score, geom="histogram", main =sprintf("%s",sample_name), xlab = "BCDS score", ylab="Counts",binwidth=0.05)
  
}

list <-lapply(sample_names, make_doublet_hist)

png("quality_control/doublet_histograms.png", width = 3000, height = 2000, units="px")
cowplot::plot_grid(plotlist = list)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


