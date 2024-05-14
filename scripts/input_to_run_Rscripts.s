#!/bin/bash
#Maria Mikedis
#maria.mikedis@cchmc.org


# link CellRanger output
CELLRANGER="210416_WIGTC-NOVASEQ1A_AH272KDSX2/cellranger"  
CELLRANGER2="210709_WIGTC-NOVASEQ1B_BHCWLNDSX2/cellranger"  

DIR="MM337_CellRanger_quant"
cd ${DIR}
mkdir filtered_feature_bc_matrix
cd filtered_feature_bc_matrix


for SAMPLE in L16_2807 L16_2808; do
	cd ${DIR}/filtered_feature_bc_matrix
	mkdir ${SAMPLE}
	ln -s ${CELLRANGER}/${SAMPLE}/filtered_feature_bc_matrix/barcodes.tsv.gz ${SAMPLE}/barcodes.tsv.gz
	ln -s ${CELLRANGER}/${SAMPLE}/filtered_feature_bc_matrix/features.tsv.gz ${SAMPLE}/features.tsv.gz	
	ln -s ${CELLRANGER}/${SAMPLE}/filtered_feature_bc_matrix/matrix.mtx.gz ${SAMPLE}/matrix.mtx.gz
done



for SAMPLE in  L16_3018 L16_3020; do
	cd ${DIR}/filtered_feature_bc_matrix
	mkdir ${SAMPLE}
	ln -s ${CELLRANGER2}/${SAMPLE}/filtered_feature_bc_matrix/barcodes.tsv.gz ${SAMPLE}/barcodes.tsv.gz
	ln -s ${CELLRANGER2}/${SAMPLE}/filtered_feature_bc_matrix/features.tsv.gz ${SAMPLE}/features.tsv.gz	
	ln -s ${CELLRANGER2}/${SAMPLE}/filtered_feature_bc_matrix/matrix.mtx.gz ${SAMPLE}/matrix.mtx.gz
done

cd ${DIR}

## quality control
bsub Rscript quality_control_pipeline.R 


## filter and integrate
bsub Rscript filter_and_integrate_data_WT.MeiocKO.R 



### identify cell types in WT samples only
bsub Rscript assign_cell_types_WT.samples.only_celltype.markers.only.R ## identify clusters based on marker gene expression; analyze cell cycle phase
bsub Rscript assign_cell_types_WT.samples.only_markers.and.cellcycle.R	# identify clusters by integrating marker gene expression and cell cycle phase

### using cell types determined using WT samples only, assign celltypes to both WT and Meioc KO samples
bsub Rscript assign_cell_types_WT.MeiocKO.R




### pseudotime analysis
bsub Rscript pseudotime_convert.from.seu_WT.MeiocKO.R
bsub Rscript pseudotime_WT.MeiocKO_graph.R

