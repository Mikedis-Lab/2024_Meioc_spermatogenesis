#Maria Mikedis
#maria.mikedis@cchmc.org



#Script assigns cell types to wild-type and Meioc knockout cells. It assigns cell types via marker expression and cell cycle designations and merges clusters as needed based on wild-type assignments worked out in assign_cell_types_WT.samples.only_celltype.markers.only.R and assign_cell_types_WT.samples.only_markers.and.cellcycle.R. 


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






###READ IN SEURAT OBJECT###
setwd(base_dir)
seu <- readRDS(file=rds_file)
DefaultAssay(seu) <- "integrated" 




###LABEL CLUSTERS
#Rename identities using new cluster IDs
names(new.cluster.ids) <- levels(seu)

seu <- RenameIdents(seu, new.cluster.ids)
seu$celltype <- Idents(seu)





#Order cell types for labeled umap and violin plots, with germ cell trajectory first, then somatic cells
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA")
                )




my_levels <- na.omit(unique(new.cluster.ids))

my_levels =  na.omit(factor(my_levels,
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
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
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", 
                 "#7f96ff", ## mut spermatocytes color
                 "#878787",  "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF","#D3D3D3"
                 # "#E38900","#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )



show_col(color_list)### to see colors




#Make new labeled UMAP

tiff("umap_labeled_integrated_original_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list, repel=FALSE) + 
  NoLegend() + 
  theme_minimal
dev.off()


tiff("umap_labeled_integrated_original_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()



#seu <- readRDS(file="integrated_data_processed_labeled_cellcycle.rds")


### with Meioc KO cells only
seu.WT = subset(seu, subset = (Genotype %in% c("WT")))
seu.WT.trim = subset(seu.WT, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  #"Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro" #,"Soma UA"
                                  )))



#### Meioc KO cells only, all clusters (including unassigned clusters)
### made my own color list from the above code; color list optimized for color distribution across clusters, without Germ UA and Soma UA
 color_list1 = c("#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", 
                 "#7f96ff", ## mut spermatocytes color
                 "#878787",  #Germ USA
                 "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF" ,"#D3D3D3" #last one is Soma UA
                 # "#E38900","#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )

tiff("umap_labeled_integrated_original_WTonly_all.clusters_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu.WT, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list1, repel=FALSE) + 
  NoLegend() + 
  scale_x_continuous(limits = c(-11, 15), breaks = c(-10, -5, 0, 5, 10, 15)) +
  theme_minimal
dev.off()


tiff("umap_labeled_integrated_original_WTonly_all.clusters_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu.WT, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list1, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()







### with Meioc KO cells only
seu.KO = subset(seu, subset = (Genotype %in% c("MeiocKO")))
seu.KO.trim = subset(seu.KO, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  #"Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro" #,"Soma UA"
                                  )))



#### Meioc KO cells only, all clusters (including unassigned clusters)
### made my own color list from the above code; color list optimized for color distribution across clusters, without Germ UA and Soma UA
 color_list1 = c("#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", 
                 "#7f96ff", ## mut spermatocytes color
                 "#878787",  #Germ USA
                 "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF" ,"#D3D3D3" #last one is Soma UA
                 # "#E38900","#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )

tiff("umap_labeled_integrated_original_MeiocKOonly_all.clusters_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu.KO, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list1, repel=FALSE) + 
  NoLegend() + 
  scale_x_continuous(limits = c(-11, 15), breaks = c(-10, -5, 0, 5, 10, 15)) +
  theme_minimal
dev.off()


tiff("umap_labeled_integrated_original_MeiocKOonly_all.clusters_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu.KO, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list1, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()


### Meioc KO cells only, clusters exluding unassinted clusters
### made my own color list from the above code; color list optimized for color distribution across clusters, without Germ UA and Soma UA
 color_list2 = c("#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8", 
                 "#7f96ff", ## mut spermatocytes color
                 #"#878787",  
                 "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF" #,"#D3D3D3"
                 # "#E38900","#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )

tiff("umap_labeled_integrated_original_MeiocKOonly_unlabeled.tiff", units="in", width=3, height=2.5, res=1200)
  DimPlot(seu.KO.trim, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list2, repel=FALSE) + 
  NoLegend() + 
  scale_x_continuous(limits = c(-11, 15), breaks = c(-10, -5, 0, 5, 10, 15)) +
  theme_minimal
dev.off()


tiff("umap_labeled_integrated_original_MeiocKOonly_unlabeled_legend.tiff", units="in", width = 8, height = 5, res=1200)
  DimPlot(seu.KO.trim, reduction = "umap", label = FALSE, label.size=4, pt.size = 0.5, cols=color_list2, repel=FALSE) + 
  #NoLegend() + 
  theme_minimal
dev.off()







###SAVE LABELED RDS FILE###
saveRDS(seu, file = "integrated_data_processed_labeled.rds")





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
saveRDS(seu, file = "integrated_data_processed_labeled_cellcycle.rds")






pdf("umap_labeled_integrated_original_cell.cycle_legend.pdf", width = 4.5, height = 2.5)
  DimPlot(seu, reduction = "umap", label = FALSE, group.by="Phase", label.size=4, pt.size = 0.5, repel=FALSE) + 
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
## WT cells only
cellcycle.WT = data.frame()
seu.WT = subset(seu, subset = (Genotype %in% "WT"))
for (group in unique(seu.WT@meta.data$celltype)) { 
    seu.subset <- subset(seu.WT, subset = (celltype %in% group))
    tmp = table(seu.subset@meta.data$Phase)
    tmp2 =data.frame(tmp)
    tmp3 = data.frame(tmp2$Freq)
    row.names(tmp3) = tmp2$Var1
    colnames(tmp3) = c(group)
    if (nrow(cellcycle.WT) == 0) {
      cellcycle.WT = tmp3
    } else {
      cellcycle.WT = cbind(cellcycle.WT, tmp3)
    }
}

cellcycle.WT = t(as.matrix(cellcycle.WT))
cellcycle.WT.ratio = as.matrix(cellcycle.WT/rowSums(cellcycle.WT)) ### doing the math with colSums didn't work properly
cellcycle.WT.percent = t(as.matrix(cellcycle.WT.ratio))*100


col.order = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  ))
cellcycle.WT = cellcycle.WT[rev(col.order),]
cellcycle.WT.percent = cellcycle.WT.percent[,col.order]

row.order = c("G1", "S", "G2M")
cellcycle.WT.percent = cellcycle.WT.percent[row.order,]


#Make new directory to hold violin plots
if (!dir.exists("cell_cycle")) {
  system(sprintf("mkdir cell_cycle"))
}

#graph of cell cycle distribution for each cell type
colors = c("#b3e2cd", "#fdcdac", "#cbd5e8")

pdf("cell_cycle/cell_cycle_ratio_per_celltype_WT.only.pdf", width = 2.5, height = 4.5)
 par(las=1)  # keep axis label horizontal
barplot(cellcycle.WT.percent,  
        legend=row.names(cellcycle.WT.percent), 
        col=colors, 
        #beside=T, 
        cex.names = 0.5,
        cex.axis = 0.5,
        horiz=TRUE,
        border="white", 
        xlab="Percent")
dev.off()

write.table(cellcycle.WT, "cell_cycle/cell_cycle_no.cells.per.cluster_WT.only.txt", sep="\t", quote=F)


### cell cycle analysis, KO cells only
cellcycle.KO = data.frame()
seu.KO = subset(seu, subset = (Genotype %in% "MeiocKO"))

for (group in unique(seu.KO@meta.data$celltype)) { 
    seu.subset <- subset(seu.KO, subset = (celltype %in% group))
    tmp = table(seu.subset@meta.data$Phase)
    tmp2 =data.frame(tmp)
    tmp3 = data.frame(tmp2$Freq)
    row.names(tmp3) = tmp2$Var1
    colnames(tmp3) = c(group)
    if (nrow(cellcycle.KO) == 0) {
      cellcycle.KO = tmp3
    } else {
      cellcycle.KO = cbind(cellcycle.KO, tmp3)
    }
}

cellcycle.KO = t(as.matrix(cellcycle.KO))
cellcycle.KO.ratio = as.matrix(cellcycle.KO/rowSums(cellcycle.KO)) ### doing the math with colSums didn't work properly
cellcycle.KO.percent = t(as.matrix(cellcycle.KO.ratio))*100


col.order = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  ))
cellcycle.KO = cellcycle.KO[rev(col.order),]
cellcycle.KO.percent = cellcycle.KO.percent[,col.order]

row.order = c("G1", "S", "G2M")
cellcycle.KO.percent = cellcycle.KO.percent[row.order,]


pdf("cell_cycle/cell_cycle_ratio_per_celltype_KO.only.pdf", width = 2.5, height = 4.5)
 par(las=1)  # keep axis label horizontal
barplot(cellcycle.KO.percent,  
        legend=row.names(cellcycle.KO.percent), 
        col=colors, 
        #beside=T, 
        cex.names = 0.5,
        cex.axis = 0.5,
        horiz=TRUE,
        border="white", 
        xlab="Percent")
dev.off()

write.table(cellcycle.KO, "cell_cycle/cell_cycle_no.cells.per.cluster_KO.only.txt", sep="\t", quote=F)





################
### WT and KO cells on one graph per cell cycle phase
################
cellcycle.WT = read.table("cell_cycle/cell_cycle_no.cells.per.cluster_WT.only.txt", sep="\t")

cellcycle.KO = read.table("cell_cycle/cell_cycle_no.cells.per.cluster_KO.only.txt", sep="\t")



#WT samples, calculate percentages per cell cycle phase per celltype
cellcycle.WT.perc = cellcycle.WT/rowSums(cellcycle.WT)*100

#KKO samples, calculate percentages per cell cycle phase per celltype
cellcycle.KO.perc = cellcycle.KO/rowSums(cellcycle.KO)*100

# data wrangling for graphing


cellcycle.WT.perc$genotype = "WT"
cellcycle.WT.perc$celltype = row.names(cellcycle.WT.perc)


cellcycle.KO.perc$genotype = "KO"
cellcycle.KO.perc$celltype = row.names(cellcycle.KO.perc)

cellcycle = rbind(cellcycle.WT.perc, cellcycle.KO.perc)
cellcycle.trim = cellcycle[cellcycle$celltype %in% c("Undiff",  "A1-4", "In/B","B G2/M","pL G1", "pL eS", "pL lS",
                                                        "L","Z", "P"),]

### order data
cellcycle.trim$celltype = factor(cellcycle.trim$celltype, levels=c("Undiff",  "A1-4", "In/B","B G2/M","pL G1", "pL eS", "pL lS",
                                                        "L","Z", "P")
                                  )

cellcycle.trim$genotype = factor(cellcycle.trim$genotype, levels = c("WT", "KO"))


###G2M graph

tiff(sprintf(paste("cell_cycle/cell.cycle_perc.cells.in.G2M.per.celltype.tif", sep="")), units="in",width = 2.5, height = 1.5, res=1200)
    print({
        ggplot(data= cellcycle.trim, aes(x=celltype, y=G2M, group=genotype)) +
              geom_line(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2)))+ ## WT = purple and KO = green
              geom_point(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2))) + 
              ylab("Percentage") + 
              ylim(c(0,100)) +
              theme_minimal + 
              theme(axis.text.x=element_text(colour="black"), 
                    axis.text.y=element_text(colour="black")) 
          })
dev.off()

###S graph

tiff(sprintf(paste("cell_cycle/cell.cycle_perc.cells.in.S.per.celltype.tif", sep="")), units="in",width = 2.25, height = 1.5, res=1200)
    print({
        ggplot(data= cellcycle.trim, aes(x=celltype, y=S, group=genotype)) +
              geom_line(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2)))+ ## WT = purple and KO = green
              geom_point(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2))) + 
              ylab("Percentage") + 
              ylim(c(0,100)) +
              theme_minimal + 
              theme(axis.text.x=element_text(colour="black"), 
                    axis.text.y=element_text(colour="black")) 
          })
dev.off()

###G1 graph


tiff(sprintf(paste("cell_cycle/cell.cycle_perc.cells.in.G1.per.celltype.tif", sep="")), units="in",width = 2.25, height = 1.5, res=400)
    print({
        ggplot(data= cellcycle.trim, aes(x=celltype, y=G1, group=genotype)) +
          geom_rect(
            aes(xmin = celltype[2],
                xmax = celltype[3],
                ymin = - Inf,
                ymax = Inf,
                fill = "lightgrey"),
                alpha = 0.5) +
          geom_rect(
            aes(xmin = celltype[4],
                xmax = celltype[5],
                ymin = - Inf,
                ymax = Inf,
                fill = "lightgrey"),
                alpha = 0.5) +
          geom_rect(
            aes(xmin = celltype[6],
                xmax = celltype[7],
                ymin = - Inf,
                ymax = Inf,
                fill = "lightgrey"),
                alpha = 0.5) +
          geom_rect(
            aes(xmin = celltype[8],
                xmax = celltype[9],
                ymin = - Inf,
                ymax = Inf,
                fill = "lightgrey"),
                alpha = 0.5) +
              geom_line(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2)))+ ## WT = purple and KO = green
              geom_point(color=c(rep("#a48ea1",nrow(cellcycle.trim)/2),rep("#b3d393",nrow(cellcycle.trim)/2))) + 
              ylab("Percentage") + 
              ylim(0,100) +
              theme_minimal + 
              theme(axis.text.x=element_text(colour="black"), 
                    axis.text.y=element_text(colour="black"),
                     legend.position="none") 
          })
dev.off()


celltype = c("Undiff",  "A1-4", "In/B","B G2/M", "pL G1", "pL eS", "pL lS", "L","Z", "P")
stats = data.frame()
i=1
while (i<=length(celltype)){
    print(celltype[i])

    ko.numbers=unname(unlist(cellcycle.KO[row.names(cellcycle.KO) %in% celltype[i],]))
    wt.perc= unname(unlist(cellcycle.WT.perc[row.names(cellcycle.WT.perc) %in% celltype[i],][,1:3]))
    wt.ratios = wt.perc/100
    

    print(ko.numbers)
    print(wt.ratios)
    p.val = chisq.test(x=ko.numbers, p=wt.ratios)$p.val

    p.adj = p.adjust(p.val, method = "bonferroni", n = length(celltype))

    if (i==1){
      row = as.data.frame(cbind(celltype[i], p.adj))
      stats = row
    } else {
      row = as.data.frame(cbind(celltype[i], p.adj))
      stats = rbind(stats, row)
    }
    i=i+1
}

colnames(stats)[1] = "celltype"

colnames(stats)[2] = "p.adj.chisq.test"

write.table(stats, "cell_cycle/cell.cycle_perc.cells.per.celltype_stats.txt", sep="\t", quote=F, row.names=F)





###MAKE VIOLIN PLOTS SHOWING EXPRESSION OF MARKER GENES IN EACH CELL TYPE###

#Make new directory to hold violin plots
if (!dir.exists("violin_plots")) {
  system(sprintf("mkdir violin_plots"))
}

if (!dir.exists("dotplots")) {
  system(sprintf("mkdir dotplots"))
}





make_violin <- function(gene) {
  png(sprintf("violin_plots/%s.png",gene), width = 1200, height = 600, units = "px")
  print({
    VlnPlot(seu, features = c(gene), pt.size=0.05, group.by = "celltype", cols=color_list)  ### colors same as umap
  })
  dev.off()
}



for (gene in marker_genes) {
  make_violin(gene)
}



###MAKE SUMMARY DOT PLOT OF MARKER GENES###
DefaultAssay(seu) <- "RNA"


## reverse cell type order so that it appears the way I want on dotplot
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                                )
                )


seu@meta.data$Genotype=   factor(seu@meta.data$Genotype, 
                                levels = rev(c("WT", "MeiocKO")))

pdf("dot_plot_integrated.pdf", width = 4.5, height = 3)
DotPlot(seu, assay ="RNA", features = marker_genes_for_pub, cols = c("#e0f3db", "#43a2ca"), dot.scale = 2, group.by = "celltype")  + 
theme_minimal + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ## rotate x axis labels 45degress
xlab("Genes") + 
ylab("") # no axis label
dev.off()







## reverse cell type order so that it appears the way I want on dotplot
seu.WT@meta.data$celltype =   factor(seu.WT@meta.data$celltype, 
                                levels = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                                )
                )

###scale=F is needed to get WT and KO dot plots to print on the same expression color scale
pdf("dot_plot_integrated_WT.only.pdf", width = 4, height = 2.25)
DotPlot(seu.WT, assay ="RNA", features = marker_genes_for_pub, cols = c("#e0f3db", "#43a2ca"), scale=F, dot.scale = 2, group.by = "celltype")  + 
theme_minimal + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ## rotate x axis labels 45degress
xlab("Genes") + 
ylab("") # no axis label
dev.off()

## reverse cell type order so that it appears the way I want on dotplot
seu.KO@meta.data$celltype =   factor(seu.KO@meta.data$celltype, 
                                levels = rev(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                                )
                )

###scale=F is needed to get WT and KO dot plots to print on the same expression color scale
pdf("dot_plot_integrated_MeiocKO.only.pdf", width = 4, height = 2.25)
DotPlot(seu.KO, assay ="RNA", features = marker_genes_for_pub, cols = c("#e0f3db", "#43a2ca"),scale=F,dot.scale = 2, group.by = "celltype")  + 
theme_minimal + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ## rotate x axis labels 45degress
xlab("Genes") + 
ylab("") # no axis label
dev.off()






## function to create dataframe with relative ratios of cell cycle phases per cell type
feature = data.frame()
#seu.subset = subset(seu, subset = (celltype %in% "pL G1"))
for (group in c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro"  )) { 
    print(group)
    seu.subset <- subset(seu, subset = (celltype %in% group))
    WT <-  subset(seu.subset, subset = (Genotype %in% c("WT")))
    KO <- subset(seu.subset, subset = (Genotype %in% c("MeiocKO")))
    p.val = wilcox.test(WT$nFeature_RNA, KO$nFeature_RNA )$p.val
    adj.p = p.adjust(p.val, method="bonferroni", length(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro")))
    tmp =data.frame(group, p.val, adj.p)
    if (nrow(feature) == 0) {
      feature = tmp
    } else {
      feature = rbind(feature, tmp)
    }
}
write.table(feature, file = "violin_plots/nFeature_RNA.per.geno_wilcoxtest.txt", append = FALSE, quote = F, sep = "\t", row.names = TRUE, col.names = TRUE)

# seu <- readRDS(file="integrated_data_processed_labeled_cellcycle.rds")



## cell type order for violing plots
seu@meta.data$celltype =   factor(seu@meta.data$celltype, 
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA",
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                )



dodge <- position_dodge(width = 1)


colors.WT.v.KO = c("#c6b8c4", "#d9e9c9")
## oder genotype
seu@meta.data$Genotype =   factor(seu@meta.data$Genotype, 
                                levels = c("WT", "MeiocKO"
                              ))




pdf(sprintf("violin_plots/nFeature_RNA_genotype.pdf"), width = 8.5, height = 2.5)
  VlnPlot(seu, features = c("nFeature_RNA"), pt.size=0, split.by = "Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
  #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
    coord_cartesian(ylim=c(0, 7800)) + 
  theme_minimal
dev.off()


pdf(sprintf("violin_plots/nCount_RNA_genotype.pdf"), width = 8.5, height = 2.5)
  VlnPlot(seu, features = c("nCount_RNA"), pt.size=0, split.by = "Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
 #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1)  + 
  theme_minimal
 # coord_cartesian(ylim=c(0, 150000)) ### zoom in axes without throwing out datapoints
dev.off()


pdf(sprintf("violin_plots/nCount_RNA_transformed_genotype.pdf"), width = 8.5, height = 2.5)
  VlnPlot(seu, features = c("nCount_RNA"), pt.size=0, split.by = "Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
  #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
   scale_y_log10() + 
  theme_minimal
 # coord_cartesian(ylim=c(0, 150000)) ### zoom in axes without throwing out datapoints
dev.off()



dodge <- position_dodge(width = 0.9)
## nFeature and nCount plots, one per genotype, GCs only

seu.GC = subset(seu, subset = (celltype %in% c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA"#,
                                  #"Sert", "Ley", "Myoid", "Endoth","Macro" #,"Soma UA"
                                  )))

## cell type order for violing plots
seu.GC@meta.data$celltype =   factor(seu.GC@meta.data$celltype, 
                                levels = c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Mut",
                                  "Germ UA"#,
                                  #"Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                )

seu.GC@meta.data$Genotype =   factor(seu.GC@meta.data$Genotype, levels=c("WT", "MeiocKO"))


pdf(sprintf("violin_plots/nFeature_RNA_GermCellsonly_byGenotype.pdf"), width = 7.5, height = 2)
  VlnPlot(seu.GC, features = c("nFeature_RNA"), pt.size=0, split.by="Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
  #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, 
              outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
    coord_cartesian(ylim=c(0, 7800)) + 
  theme_minimal
dev.off()


pdf(sprintf("violin_plots/nCount_RNA_GermCellsonly_byGenotype.pdf"), width = 7.5, height = 2)
  VlnPlot(seu.GC, features = c("nCount_RNA"), pt.size=0, split.by="Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
 #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, 
                outlier.size = -1, lwd=0.2)  + 
     scale_y_log10() + 
  theme_minimal #+
 #coord_cartesian(ylim=c(0, 80000)) ### zoom in axes without throwing out datapoints
dev.off()






## nFeature and nCount plots, one per genotype, GCs only

seu.Soma = subset(seu, subset = (celltype %in% c(#"Undiff",  "A1-4",
                                  # "In/B","B G2/M",
                                  # "pL G1", "pL eS", "pL lS",
                                  # "L","Z", "P", 
                                  # "Mut",
                                  # "Germ UA"#,
                                  "Sert", "Ley", "Myoid", "Endoth","Macro" ,"Soma UA"
                                  )))

## cell type order for violing plots
seu.Soma@meta.data$celltype =   factor(seu.GC@meta.data$celltype, 
                                levels = c(#"Undiff",  "A1-4",
                                  # "In/B","B G2/M",
                                  # "pL G1", "pL eS", "pL lS",
                                  # "L","Z", "P", 
                                  # "Mut",
                                  # "Germ UA"#,
                                  "Sert", "Ley", "Myoid", "Endoth","Macro",  "Soma UA"
                                  )
                )

seu.Soma@meta.data$Genotype =   factor(seu.Soma@meta.data$Genotype, levels=c("WT", "MeiocKO"))


pdf(sprintf("violin_plots/nFeature_RNA_Somaonly_byGenotype.pdf"), width = 4.5, height = 2)
  VlnPlot(seu.Soma, features = c("nFeature_RNA"), pt.size=0, split.by="Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
  #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, 
              outlier.size = -1, lwd=0.2) + ### lwd to adjust line thickness
    coord_cartesian(ylim=c(0, 7800)) + 
  theme_minimal
dev.off()


pdf(sprintf("violin_plots/nCount_RNA_SomaCellsonly_byGenotype.pdf"), width = 4.5, height = 2)
  VlnPlot(seu.Soma, features = c("nCount_RNA"), pt.size=0, split.by="Genotype", group.by = "celltype",  cols=colors.WT.v.KO) + 
 #NoLegend() + 
  geom_boxplot(width=0.1, position = dodge, 
                outlier.size = -1, lwd=0.2)  + 
     scale_y_log10() + 
  theme_minimal #+
 #coord_cartesian(ylim=c(0, 80000)) ### zoom in axes without throwing out datapoints
dev.off()








## function to create dataframe with relative ratios of cell cycle phases per cell type
count = data.frame()
#seu.subset = subset(seu, subset = (celltype %in% "pL G1"))
for (group in c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Sert", "Ley", "Myoid", "Endoth","Macro")) { 
    print(group)
    seu.subset <- subset(seu, subset = (celltype %in% group))
    WT <-  subset(seu.subset, subset = (Genotype %in% c("WT")))
    KO <- subset(seu.subset, subset = (Genotype %in% c("MeiocKO")))
    p.val = wilcox.test(WT$nCount_RNA, KO$nCount_RNA )$p.val
    adj.p = p.adjust(p.val, method="bonferroni", length(c("Undiff",  "A1-4",
                                  "In/B","B G2/M",
                                  "pL G1", "pL eS", "pL lS",
                                  "L","Z", "P", 
                                  "Sert", "Ley", "Myoid", "Endoth","Macro")))
    tmp =data.frame(group, p.val, adj.p)
    if (nrow(count) == 0) {
      count = tmp
    } else {
      count = rbind(count, tmp)
    }
}
write.table(count, file = "violin_plots/nCount_RNA.per.geno_wilcoxtest.txt", append = FALSE, quote = F, sep = "\t", row.names = TRUE, col.names = TRUE)




