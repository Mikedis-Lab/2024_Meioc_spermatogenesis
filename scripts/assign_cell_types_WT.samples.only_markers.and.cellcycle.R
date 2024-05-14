#Maria Mikedis
#maria.mikedis@cchmc.org


#Script uses the wild-type samples only. It assigns cell types via marker expression and cell cycle designations and merges clusters as needed. 


library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)
library(hash)
library(cowplot)
library(scales)


###SET VARIABLES###
base_dir <- "MM337_CellRanger_quant/filtered_feature_bc_matrix/integrated_results_WT.MeiocKO.only"  #folder that contains results from filter_and_integrate_data.R
rds_file <- "integrated_data_processed.rds"   #RDS file created by filter_and_integrate_data.R


###READ IN SEURAT OBJECT###
setwd(base_dir)
seu <- readRDS(file=rds_file)
DefaultAssay(seu) <- "integrated" 



if (!dir.exists("WT_samples_only")) {
  system(sprintf("mkdir WT_samples_only"))
}

if (!dir.exists("WT_samples_only/celltype.by.markers.and.cellcycle")) {
  system(sprintf("mkdir WT_samples_only/celltype.by.markers.and.cellcycle"))
}


###SET MARKER GENES###        
#List of genes used in feature plots from filter_and_integrate_data.R
marker_genes <- c("Cyp11a1", "Dlk1",  ## fetal Leydig
                  "Cd14", "Adgre1", "Itgam", #macrophages
                  "Acta2", "Myh11", #peritubular myoid cells
                  "Tm4sf1", #vascular endothelium
                  "Sox9", "Cldn11" , # Sertoli cells
                  "Id4", "Gfra1", "Etv5", #SSCs
                  "Zbtb16", "Sall4", #undiff gonia
                  "Kit", "Stra8", "Ccnd2", #differentiating type A spermatogonia
                  "Kit", #diff gonia Intermiate to type B
                  "Cenpa", #diff gonia, G1 and early S pL https://www.nature.com/articles/s41422-018-0074-y
                  "Stra8", "Nacad", "Myl7",#mid and late S pL https://www.nature.com/articles/s41422-018-0074-y
                  "Gm960", #Lept https://www.nature.com/articles/s41422-018-0074-y
                  "Meiob", #Lept and Zygo, a little lower in early Pachy https://www.nature.com/articles/s41422-018-0074-y
                  "Meioc", "Gm4969", 
                  "Piwil1", "Piwil2", "Hspa2") ## late spermatocytes

marker_genes_for_pub <- c("Dazl", "Ddx4",
                  "Gfra1", "Id4", #SSCs
                  "Zbtb16", #undiff gonia
                  "Kit", #diff gonia
                  "Stra8", # pL
                  "Ccnb3", "Meiob", #Lept and Zygo, a little lower in early Pachy https://www.nature.com/articles/s41422-018-0074-y
                  "Piwil1", ## pachytene
                  "Cldn11","Sox9",  # Sertoli cells
                  "Cyp11a1", "Dlk1",  ## fetal Leydig
                  "Acta2", "Myh11", #peritubular myoid cells
                  "Icam2", "Sox17", #vascular endothelium
                  "Adgre1", "Cd14") #macrophages
        



### theme for manuscript figures
theme_minimal <- theme(text=element_text(family="Helvetica",face="plain",size=6),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=6),
                  axis.line.y=element_blank(),
                  axis.line.x=element_line(size=0.25),
                  axis.ticks=element_line(size=0.25),
                  legend.title=element_blank(),
                  legend.background=element_rect(fill="transparent"),
                  legend.text=element_text(size=6), plot.title=element_text(size=6, face="plain"),
                  panel.background = element_blank(),
                  legend.key=element_blank())


###ASSIGN CELL TYPES FOR EACH CLUSTER###
#Use files created by filter_and_integrate_data.R:  
  #1. umap_integrated.png  (shows clusters)
  #2. Plots in feature_plots/ (shows expression of marker genes)

# cell types based on markers and cell cycle analysis from seurat


## cluster ids based on marker expression and cell cycle designation
new.cluster.ids <- c( "Myoid", "pL G1", "pL lS", "L", "B G2/M", #0
						"In/B", "Z", "Z", "In/B", "Ley", #5
						"L", "Germ UA", "Undiff", "pL eS", "A1-4", #10
						"In/B", "pL lS", "Sert", "Sert", "Undiff", #15
						"P", "A1-4", "A1-4", "Myoid", "P", #20
						"Mut", "A1-4", "In/B", "Sert", "Ley", #25
						"Macro", "Endoth", "Endoth", "Macro", "Soma UA", #30
						"Sert", "Myoid") #35








###LABEL CLUSTERS
#Rename identities using new cluster IDs
names(new.cluster.ids) <- levels(seu)

seu <- RenameIdents(seu, new.cluster.ids)
seu$celltype <- Idents(seu)

## select WT cells only
seu.WT = subset(seu, subset = (Genotype %in% "WT"))

seu=seu.WT

### Drop clusters with less than 10 cells --> WT cells have 4 cells in the mutant cluster; this filter let's me drop the mutant cluster until later. 

seu2 = subset(seu, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")))



seu=seu2


#Order cell types for labeled umap and violin plots, with germ cell trajectory first, then somatic cells
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")
                )




my_levels <- na.omit(unique(new.cluster.ids))

my_levels =  na.omit(factor(my_levels,
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")
                                  ))


Idents(seu) <- factor(Idents(seu), levels= sort(my_levels))




###MAKE INTEGRATED UMAP WITH NEW CLUSTER LABELS####

#Generate list of colors for new integrated plot  
#Replace "Unassigned" with gray (#D3D3D3)





color_list <- hue_pal()(length(my_levels))
color_list = sample(color_list) # randomize colors

#color_list <- hue_pal()(length(my_levels))
for (i in 1:length(color_list)) {
  if (my_levels[i]=="Soma UA"  ) {
    color_list[i]="#D3D3D3"
  }
  if (my_levels[i]=="Germ UA"  ) {
    color_list[i]="#878787"
  }
}



### made my own color list from the above code; color list optimized for color distribution across clusters
 color_list = c("#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", "#878787",  "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF","#D3D3D3"
                 #"#E38900", "#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )


show_col(color_list)### to see colors



#Make new labeled UMAP

tiff("WT_samples_only/celltype.by.markers.and.cellcycle/umap_labeled_integrated_original_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list, repel=FALSE) + 
  NoLegend() + 
  theme_minimal
dev.off()


tiff("WT_samples_only/celltype.by.markers.and.cellcycle/umap_labeled_integrated_original_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()





###SAVE LABELED RDS FILE###
saveRDS(seu, file = "WT_samples_only/celltype.by.markers.and.cellcycle/integrated_data_processed_labeled.rds")





### Cell cycle analysis designate cell type subpopulations
# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle (raw data file version). Read it in with:

DefaultAssay(seu) <- "RNA"

cell_cycle_genes <- read.csv("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 


library('biomaRt')

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
#listFilters(ensembl) # to see list of possible filters
mart = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
      filters = 'ensembl_gene_id', 
      values = cell_cycle_genes$geneID, 
      mart = ensembl, 
      useCache = FALSE) # useCache setting needed to avoid an error

cell_cycle_genes2 = merge(cell_cycle_genes,mart,by.x="geneID",by.y="ensembl_gene_id")


# Acquire the S phase genes
s_genes <- cell_cycle_genes2 %>%
        dplyr::filter(phase == "S") %>%
        pull("external_gene_name")
        
# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_genes2 %>%
        dplyr::filter(phase == "G2/M") %>%
        pull("external_gene_name")

# Perform cell cycle scoring
seu <- CellCycleScoring(seu,
                        g2m.features = g2m_genes,
                        s.features = s_genes)



# order cell cycle features
seu$Phase = factor(seu$Phase, levels= c("G1", "S", "G2M"))

###SAVE LABELED RDS FILE###
saveRDS(seu, file = "WT_samples_only/celltype.by.markers.and.cellcycle/integrated_data_processed_labeled_cellcycle.rds")


#seu <- readRDS(file="WT_samples_only/celltype.by.markers.and.cellcycle/integrated_data_processed_labeled_cellcycle.rds")



pdf("umap_labeled_integrated_original_cell.cycle_legend.pdf", width = 4.5, height = 2.5)
  DimPlot(seu, reduction = "umap", label = FALSE, group.by="Phase", label.size=4, pt.size = 0.5, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()

png("umap_labeled_integrated_original_cell.cycle_genotype_legend.png", width = 1200, height = 600, units = "px")
  DimPlot(seu, reduction = "umap", label = FALSE, group.by="Phase", , split.by="Genotype", label.size=4, pt.size = 0.5, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()

pdf("umap_labeled_integrated_cell.cycle_unlabeled.pdf", width = 3.5, height = 2.5)
  DimPlot(seu, reduction = "umap", label = FALSE, group.by="Phase", label.size=4, pt.size = 0.5, repel=FALSE) +   
  NoLegend() + 
  ggtitle("") + 
  theme_minimal
dev.off()






## function to create dataframe with relative ratios of cell cycle phases per cell type
cellcycle = data.frame()
#seu.WT = subset(seu, subset = (Genotype %in% "WT"))
for (group in unique(seu@meta.data$celltype)) { 
    seu.subset <- subset(seu, subset = (celltype %in% group))
    tmp = table(seu.subset@meta.data$Phase)
    tmp2 =data.frame(tmp)
    tmp3 = data.frame(tmp2$Freq)
    row.names(tmp3) = tmp2$Var1
    colnames(tmp3) = c(group)
    if (nrow(cellcycle) == 0) {
      cellcycle = tmp3
    } else {
      cellcycle = cbind(cellcycle, tmp3)
    }
}

cellcycle = t(as.matrix(cellcycle))
cellcycle.ratio = as.matrix(cellcycle/rowSums(cellcycle)) ### doing the math with colSums didn't work properly
cellcycle.percent = t(as.matrix(cellcycle.ratio))*100


col.order = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  ))
cellcycle.percent = cellcycle.percent[,col.order]

row.order = c("G1", "S", "G2M")
cellcycle.percent = cellcycle.percent[row.order,]


#Make new directory to hold violin plots
if (!dir.exists("WT_samples_only/celltype.by.markers.and.cellcycle/cell_cycle")) {
  system(sprintf("mkdir WT_samples_only/celltype.by.markers.and.cellcycle/cell_cycle"))
}



#graph of cell cycle distribution for each cell type
#colors = c("#b3e2cd", "#fdcdac", "#cbd5e8")
colors = c("#7ccda9","#fba061","#96aad1")

pdf("WT_samples_only/celltype.by.markers.and.cellcycle/cell_cycle/cell_cycle_ratio_per_celltype_final.pdf", width = 1.75, height = 4.5)
 par(las=1)  # keep axis label horizontal
barplot(cellcycle.percent,  
        legend=row.names(cellcycle.percent), 
        col=colors, 
        #beside=T, 
        cex.names = 0.5,
        cex.axis = 0.5,
        horiz=TRUE,
        border="white", 
        xlab="Percent")
dev.off()




###MAKE VIOLIN PLOTS SHOWING EXPRESSION OF MARKER GENES IN EACH CELL TYPE###

#Make new directory to hold violin plots
if (!dir.exists("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots")) {
  system(sprintf("mkdir WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots"))
}

if (!dir.exists("WT_samples_only/celltype.by.markers.and.cellcycle/dotplots")) {
  system(sprintf("mkdir WT_samples_only/celltype.by.markers.and.cellcycle/dotplots"))
}





make_violin <- function(gene) {
  png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/%s.png",gene), width = 1200, height = 600, units = "px")
  print({
    VlnPlot(seu, features = c(gene), pt.size=0.05, group.by = "celltype", cols=color_list)  ### colors same as umap
  })
  dev.off()
}



for (gene in marker_genes) {
  make_violin(gene)
}

  png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/Ccna2_genotype.png",gene), width = 2400, height = 600, units = "px")
    VlnPlot(seu, features = c("Ccna2"), pt.size=0.05, group.by = "celltype", split.by="Genotype")
  dev.off()

    png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/Ccna2.png",gene), width = 1200, height = 600, units = "px")
    VlnPlot(seu, features = c("Ccna2"), pt.size=0.05, group.by = "celltype")
  dev.off()

      png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/Dazl.png",gene), width = 1200, height = 600, units = "px")
    VlnPlot(seu, features = c("Dazl"), pt.size=0.05, group.by = "celltype")
  dev.off()


    png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/Ddx4.png",gene), width = 1200, height = 600, units = "px")
    VlnPlot(seu, features = c("Ddx4"), pt.size=0.05, group.by = "celltype")
  dev.off()

      png(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/Mki67.png",gene), width = 1200, height = 600, units = "px")
    VlnPlot(seu, features = c("Mki67"), pt.size=0.05, group.by = "celltype")
  dev.off()


###MAKE SUMMARY DOT PLOT OF MARKER GENES###
DefaultAssay(seu) <- "RNA"


## reverse cell type order so that it appears the way I want on dotplot
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")
                                ))
                


tiff("WT_samples_only/celltype.by.markers.and.cellcycle/dotplots/dot_plot_integrated.tiff", units="in", width=4, height=2.5, res=1200)
DotPlot(seu, assay ="RNA", features = marker_genes_for_pub, cols = c("#e0f3db", "#43a2ca"), dot.scale = 2, group.by = "celltype")  + 
theme_minimal + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ## rotate x axis labels 45degress
xlab("Genes") + 
ylab("") # no axis label
dev.off()




### dotplot with Soma UA and Germ UA
seu.trim = subset(seu, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  #"Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro"#,  "Soma UA"
                                  )))



tiff("WT_samples_only/celltype.by.markers.and.cellcycle/dotplots/dot_plot_integrated_no.unassigned.tiff", units="in", width=4, height=2.5, res=1200)
DotPlot(seu.trim, assay ="RNA", features = marker_genes_for_pub, cols = c("#e0f3db", "#43a2ca"), dot.scale = 2, group.by = "celltype")  + 
theme_minimal + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ## rotate x axis labels 45degress
xlab("Genes") + 
ylab("") # no axis label
dev.off()

### colorlist without unassigned populations
 color_list2 = c("#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", #"#878787",  
                 "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF"#,"#D3D3D3"
                 #"#E38900", "#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )


#### UMAP without unassigned populations
tiff("WT_samples_only/celltype.by.markers.and.cellcycle/umap_labeled_integrated_original_no.unassigned_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu.trim, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list2, repel=FALSE) + 
  NoLegend() + 
  theme_minimal
dev.off()


tiff("WT_samples_only/celltype.by.markers.and.cellcycle/umap_labeled_integrated_original_no.unassigned_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu.trim, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list2, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()



dodge <- position_dodge(width = 1)



## reorder cell types so that it appears in correct order on violin plot
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  #"Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")
                                )


pdf(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/nFeature_RNA.pdf"), width = 6.5, height = 1.8)
  VlnPlot(seu, features = c("nFeature_RNA"), pt.size=0, group.by = "celltype",  cols=color_list) + 
  NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
    #coord_cartesian(ylim=c(0, 7800)) + 
  theme_minimal
dev.off()


pdf(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/nCount_RNA.pdf"), width = 6.5, height = 1.8)
  VlnPlot(seu, features = c("nCount_RNA"), pt.size=0, group.by = "celltype",  cols=color_list) + 
  NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1)  + 
  theme_minimal
 # coord_cartesian(ylim=c(0, 150000)) ### zoom in axes without throwing out datapoints
dev.off()


pdf(sprintf("WT_samples_only/celltype.by.markers.and.cellcycle/violin_plots/nCount_RNA_transformed.pdf"), width = 6.5, height = 1.8)
  VlnPlot(seu, features = c("nCount_RNA"), pt.size=0, group.by = "celltype",  cols=color_list) + 
  NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
   scale_y_log10() + 
  theme_minimal
 # coord_cartesian(ylim=c(0, 150000)) ### zoom in axes without throwing out datapoints
dev.off()


