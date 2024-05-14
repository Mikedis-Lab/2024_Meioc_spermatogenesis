#Maria Mikedis
#maria.mikedis@cchmc.org


#This script will:  
#1) Filter cells based on criteria from Quality Control Pipeline 
#2) Integrate data between samples 
#3) Add in your metadata to the Seurat object (minimum batch and genotype)
#4) Make a UMAP
#5) Generate feature plots to help assign cell types

library(BUSpaRse)
library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(DropletUtils)
library(Matrix)
theme_set(theme_bw())
library(reticulate)
library(dplyr)
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(hash)
library(reshape)


###SET FILENAMES###
base_dir <- "MM337_CellRanger_quant/filtered_feature_bc_matrix/"   #directory with files matrix.mx.gz, features.tsv.gz, and barcodes.tsv.gz, one per sample
metadata <- "genome_files/Mus_musculus_Gencode.vM23.basic/gencode.vM23.basic_gcna_txo_annotable.txt" #tab-separated gene metadata file containing, at minimum, columns "gene_name" and "gene_type" (i.e "protein_coding")
setwd(base_dir)

if (!dir.exists("integrated_results_WT.MeiocKO.only")) {
  system(sprintf("mkdir integrated_results_WT.MeiocKO.only"))
}


###SET QC THRESHOLDS###

min_features <- 1000   ### 500 is too low
max_mito <- 10
max_doublet <- 0.4 ## strict doublet filter


###CHOOSE SAMPLES###   these should be the names of samples you want to include from the bus_output folder (excluding any outliers!) 
samples_list <- c("L16_2807", "L16_2808", "L16_3018", "L16_3020")  

###READ IN GENE METADATA###
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


###GENERATE FILTERED SEURAT OBJECT FOR EACH SAMPLE###
object_list <- c()

for (sample in samples_list) {
  
    #Read in file
  #directory <- sprintf("bus_output/%s/genecount", sample)
  res_mat <- Read10X(data.dir = paste(base_dir, sample, sep=""))
  
  #Filter to protein-coding genes only
  res_mat <- res_mat[row.names(res_mat) %in% genes_of_interest,]
  
  #Make object
  object <- CreateSeuratObject(counts=res_mat, min.cells = 3, min.features=min_features, project=sample) 
  
  #Find fractions mitochondrial and ribosomal
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^Rp[ls]")
  
  #Add doublet scores
  sce <- as.SingleCellExperiment(object)
  sce <- cxds(sce,retRes = TRUE)  #coexpression
  sce <- bcds(sce,retRes = TRUE,verb=TRUE) #artificial doublets
  sce <- cxds_bcds_hybrid(sce)
  
  object@meta.data$hybrid_score <- sce$hybrid_score
  object@meta.data$bcds_score <- sce$bcds_score
  object@meta.data$cxds_score <- sce$cxds_score

  #Filter and normalize object
  object <- subset(object, subset = percent.mt < max_mito & bcds_score <max_doublet) 
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object<- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  
  #Add to list
  object_list <- c(object_list,object)
  
}

###INTEGRATE DATA###
seu.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, reference = c(1)) #takes around 15 min for 4 samples
seu <- IntegrateData(anchorset = seu.anchors, dims = 1:30)
DefaultAssay(seu) <- "integrated"

###READ IN SAMPLE METADATA### 
meta <- read.table(file="metadata.txt", header=T)

hash.batch <- hash(meta$Sample_ID, meta$Batch)
hash.geno <- hash(meta$Sample_ID, meta$Genotype)

batches <- c()
genotypes <- c()

for (i in rownames(seu@meta.data)) {
  batch <- hash.batch[[seu@meta.data[i,"orig.ident"]]]
  batches <- c(batches, as.character(batch))
  genotype <- hash.geno[[seu@meta.data[i,"orig.ident"]]]
  genotypes <- c(genotypes, as.character(genotype))
}

seu <- AddMetaData(seu, batches, col.name = "Batch")
seu <- AddMetaData(seu, genotypes, col.name = "Genotype")



#check what filtered data looks like
png("integrated_results_WT.MeiocKO.only/features_versus_mito_filtered.png", width = 1480, height = 740, units = "px")
FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=0.5) + 
  geom_hline(yintercept=15, linetype="dashed", color = "blue", size=0.5) + 
  geom_vline(xintercept=1000, linetype="dashed", color = "black", size=0.5) + 
  geom_vline(xintercept=500, linetype="dashed", color = "blue", size=0.5) 
dev.off()




# Basic function to convert human to mouse gene names
# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
# site accessed 5/21/21
ConvertHumanGeneListToMM <- function(x){
  
  # Load human ensembl attributes
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Load mouse ensembl attributes
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  
  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])
  
  # # Print the first 6 genes found to the screen
  # print(head(mouse.gene.list))
  return(mouse.gene.list)
}


###ADD CELL CYCLE INFO###    https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html
#s.genes <- cc.genes$s.genes
#s.genes.mouse = ConvertHumanGeneListToMM (s.genes) # original 43 human genes, converts to 40 mouse genes (i.e., some genes missing)
##g2m.genes <- cc.genes$g2m.genes
#g2m.genes.mouse = ConvertHumanGeneListToMM (g2m.genes)	# originally 54 human genes, converts to 49 mouse genes (i.e., some genes missing)
#seu <- CellCycleScoring(seu, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse, set.ident = TRUE)
#seu$CC.Difference <- seu$S.Score - seu$G2M.Score

###SAVE OBJECT###

saveRDS(seu, file = "integrated_results_WT.MeiocKO.only/integrated_data.rds")

###SCALE DATA AND RUN PCA###
set.seed(8)

all.genes <- rownames(seu)

if (length(unique(seu@meta.data$Batch))>1) {
  seu <- ScaleData(seu, verbose=FALSE,vars.to.regress = c("percent.mt", "Batch")) #, "CC.Difference")) ### tutorial also includes "nCount_RNA" https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html#scaling-the-data-and-removing-unwanted-sources-of-variation
} else {
  seu <- ScaleData(seu, verbose=FALSE,vars.to.regress = c("percent.mt")) #, "CC.Difference"))
}

  
set.seed(8)
seu <- RunPCA(seu, verbose = FALSE, npcs = 100)

#Find clusters and run UMAP
set.seed(8)
seu <- FindNeighbors(seu, reduction = "pca", dims=1:30,  nn.eps = 0.5) #lower = more exact, slower:  default = 0 #original used dims=1:75
set.seed(8)
seu <- FindClusters(seu, resolution = 1.5, graph.name = 'integrated_snn')
set.seed(8)
seu <- RunUMAP(seu, min.dist = 0.3, dims=1:30, metric="euclidean", n.neighbors=20L)   #original used dims=1:75


###MAKE INTEGRATED UMAP###
png("integrated_results_WT.MeiocKO.only/umap_integrated.png", width = 800, height = 600, units = "px")
DimPlot(seu, reduction = "umap", label=TRUE, repel=TRUE)
dev.off()


# Plot the elbow plot
png("integrated_results_WT.MeiocKO.only/ElbowPlot.png", width = 1480, height = 740, units = "px")
ElbowPlot(object = seu, 
          ndims = 40)
dev.off()


###SAVE FINAL FILE###
saveRDS(seu, file = "integrated_results_WT.MeiocKO.only/integrated_data_processed.rds")




#
###MAKE FEATURE PLOTS FOR CELL TYPE ASSIGNMENT###
DefaultAssay(seu) <- "RNA"

if (!dir.exists("integrated_results_WT.MeiocKO.only/feature_plots")) {
  system(sprintf("mkdir integrated_results_WT.MeiocKO.only/feature_plots"))
}
#for more germ cell markers: https://www.cell.com/cell-reports/pdf/S2211-1247(18)31602-4.pdf
png("integrated_results_WT.MeiocKO.only/feature_plots/Sertoli_cells.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Sox9", "Amh", "Cldn11"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Sertoli_cells_Hammoudpaper.png", width = 1500, height = 1000, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Clu", "Amhr2", "Sox9", "Ctsl", "Rhox8"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Innate_lymph_cells.png", width = 1500, height = 1000, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Id2","Il7r","Rora", "Thy1", "Ccl5", "Cd52"))
dev.off()



png("integrated_results_WT.MeiocKO.only/feature_plots/Peritubular_myoid_cells.png", width = 2000, height = 500, units = "px")
FeaturePlot(seu, ncol = 4, features = c("Acta2","Myh11","Myl6", "Pdgfrb"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Leydig_cells_neonatal.png", width = 500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Dlk1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Leydig_cells_adult.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Fabp3"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Leydig_cells_adult_Hammoudpaper.png", width = 2000, height = 500, units = "px")
FeaturePlot(seu, ncol = 4, features = c("Cyp17a1", "Cyp11a1", "Star", "Hsd3b1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Vascular_endothelium.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Tm4sf1","Vegf","Pecam"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Vascular_endothelium_hammoudpaper.png", width = 2000, height = 500, units = "px")
FeaturePlot(seu, ncol = 4, features = c("Vwf","Tie1","Tek", "Ly6c1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Macrophages1.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Cd14", "Adgre1", "Itgam"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Macrophages2.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Aifm1", "Cx3cr1", "Fcgr1", "H2"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Macrophages_hammoudpaper.png", width = 2000, height = 500, units = "px")
FeaturePlot(seu, ncol = 4, features = c("Apoe", "Dab2", "Cd74", "Adgre1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatogonial_stem_cells.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Id4","Gfra1","Etv5"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatogonia.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Dmrt1","Zbtb16","Kit"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatogonia_diff_A1.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Ccnd2","Stra8","Kit"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatocytes_early.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Stra8","Meioc","Gm4969"))
dev.off()

# for markers of pL vs. L vs. Z
png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatocytes_early2.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Nacad","Myl7","Gm960"))
dev.off()

# for markers of pL vs. L vs. Z
png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatocytes_early3.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Meiob","Psma8"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Spermatocytes_late.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Piwil1","Piwil2","Hspa2"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Round_spermatids.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Tex21","Zbtb16","Kit"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Elongating_spermatids.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Tnp1","Zbtb16","Kit"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/G2M_Cenpa.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, features = c("Cenpa"), split.by="Genotype")
dev.off()


png("integrated_results_WT.MeiocKO.only/feature_plots/quiescence.png", width = 500, height = 500, units = "px")
FeaturePlot(seu, features = c("Mki67"))
dev.off()


png("integrated_results_WT.MeiocKO.only/feature_plots/S_1.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Ccne2","Mcm2","Pold3"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/S_2.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Rad51","Pcna","Rad51ap1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/G2M_1.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Cdk1","Cenpa","Bub1"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/G2M_2.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Ccnb2","G2e3","Cenpe"))
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/G2M_3.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Mki67","Aurka","Ube2c"))
dev.off()


png("integrated_results_WT.MeiocKO.only/feature_plots/Germ_cells.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("Dazl","Ddx4"))
dev.off()




png("integrated_results_WT.MeiocKO.only/feature_plots/Genotype.png", width = 1500, height = 500, units = "px")
DimPlot(seu, 
        label = TRUE, 
        split.by = "Genotype")  + NoLegend()
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Batch.png", width = 1000, height = 500, units = "px")
DimPlot(seu, 
        label = TRUE, 
        split.by = "Batch")  + NoLegend()
dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Sample.png", width = 3000, height = 500, units = "px")
DimPlot(seu, 
        label = TRUE, 
        split.by = "orig.ident")  + NoLegend()
dev.off()


dev.off()

png("integrated_results_WT.MeiocKO.only/feature_plots/Doublet_score_per_sample.png", width = 3000, height = 500, units = "px")
FeaturePlot(seu, features = c("bcds_score"), split.by="orig.ident", keep.scale = "feature")
dev.off()

