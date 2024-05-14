#Maria Mikedis
#maria.mikedis@cchmc.org

#This script will generates pseudotime graphs from previously generated monocle object. 

library(ggplot2)
library(data.table)
#library(Seurat)
#library(SeuratWrappers)
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
base_dir <- "MM337_CellRanger_quant/filtered_feature_bc_matrix/integrated_results_WT.MeiocKO.only/pseudotime_monocle_GC.only"
monocle_obj =  "integrated_data_processed_labeled_monocle_WT.MeiocKO.only.rds"


####Read in integrated Seurat object###
setwd(base_dir)
cds <- readRDS(file=monocle_obj)

## order genotypes
## reverse cell type order so that it appears the way I want on dotplot
pData(cds)$Genotype =   factor(pData(cds)$Genotype, levels=c("WT", "MeiocKO"))

#### graph WT and KO cells on separate plots
cdsWT = cds[,rownames(pData(cds)[pData(cds)$Genotype == "WT",])]
cdsKO<-cds[,rownames(pData(cds)[pData(cds)$Genotype == "MeiocKO",])]

tiff("pseudotime_WT.MeiocKO_genotype_notrajectory.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cds,color_cells_by = "Genotype",cell_size = 1, 
			group_cells_by = c("cluster", "partition"),
			show_trajectory_graph = FALSE,
	       	label_branch_points = FALSE,
	       	label_roots = FALSE,
	      	label_cell_groups=TRUE,
	      	group_label_size = 4,
	      	labels_per_group = 0,
           	label_leaves=FALSE) + 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) + 
   scale_color_manual(values = c("#a48ea1", "#b3d393"))
dev.off()



tiff("pseudotime_WT.MeiocKO_genotype_notrajectory_WTonly.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cdsWT,color_cells_by = "Genotype",cell_size = 1, 
      group_cells_by = c("cluster", "partition"),
      show_trajectory_graph = FALSE,
          label_branch_points = FALSE,
          label_roots = FALSE,
          label_cell_groups=TRUE,
          group_label_size = 4,
          labels_per_group = 0,
            label_leaves=FALSE) + 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) + 
   scale_color_manual(values = c("#a48ea1"))
dev.off()

tiff("pseudotime_WT.MeiocKO_genotype_notrajectory_KOonly.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cdsKO,color_cells_by = "Genotype",cell_size = 1, 
      group_cells_by = c("cluster", "partition"),
      show_trajectory_graph = FALSE,
          label_branch_points = FALSE,
          label_roots = FALSE,
          label_cell_groups=TRUE,
          group_label_size = 4,
          labels_per_group = 0,
            label_leaves=FALSE) + 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) + 
   scale_color_manual(values = c("#b3d393" ))
dev.off()



### cluster list used in seurat based UMAPS
 color_list = c(#"#E38900", 
  "#00AEF9", "#00BF77",
                   "#A58AFF", "#B0A100", "#7698FF",  "#FF66A8", 
                 "#00BC56","#EF7F46",  "#00BFC4", "#DF70F8",
                 "#7F96FF", ### <--- color for mutant cluster
                 "#878787",  "#F8766D","#7CAE00","#FB61D7",
                  "#D59100", "#C77CFF","#D3D3D3"
                 #, "#53B400",  "#C49A00", "#FF62C1",  "#00C094", "#F166E9", "#D3D3D3" , "#7CAE00", "#00B6EB", "#D3D3D3"
                 )

tiff("pseudotime_WT.MeiocKO_celltype.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cds,color_cells_by = "celltype",cell_size = 1,
			group_cells_by = c("cluster", "partition"),
			show_trajectory_graph = T,
	       	label_branch_points = FALSE,
	       	label_roots = FALSE,
	      	label_cell_groups=TRUE,
	      	group_label_size = 2,
	      	labels_per_group = 3,
           	label_leaves=FALSE)+ 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) + 
   scale_color_manual(values = color_list)
dev.off()

tiff("pseudotime_WT.MeiocKO_celltype_unlabeled.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cds,color_cells_by = "celltype",cell_size = 1,
      group_cells_by = c("cluster", "partition"),
      show_trajectory_graph = T,
          label_branch_points = FALSE,
          label_roots = FALSE,
          label_cell_groups=TRUE,
          group_label_size = 2,
          labels_per_group = 0,
            label_leaves=FALSE)+ 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) + 
   scale_color_manual(values = color_list)
dev.off()


tiff("pseudotime_WT.MeiocKO_cellcycle_notrajectory.tif",units="in", width = 4, height = 4, res=1200)
plot_cells(cds,color_cells_by = "Phase",cell_size = 1, 
			group_cells_by = c("cluster", "partition"),
			show_trajectory_graph = FALSE,
	       	label_branch_points = FALSE,
	       	label_roots = FALSE,
	      	label_cell_groups=TRUE,
	      	group_label_size = 4,
	      	labels_per_group = 0,
           	label_leaves=FALSE)+ 
  theme(axis.line.y=element_blank()) + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-10, 12)) #+ 
   #scale_color_manual(values = c("#b3d393" ))
dev.off()




