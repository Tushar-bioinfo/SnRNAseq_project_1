This project involves the integrated analysis of single-nucleus RNA-seq data from wild-type (WT) and familial Alzheimer’s disease (FAD) brain samples to investigate differences in cell populations and gene expression patterns. Raw 10X Genomics data for both conditions were loaded and processed using the Seurat package in R. Initial steps included quality control based on mitochondrial content, followed by scaling and metadata annotation.

The WT and FAD datasets were merged and split for normalization, variable feature selection, and dimensionality reduction via PCA. Integration was performed using reciprocal PCA (RPCA) to correct batch effects and produce a unified dataset. This integrated object will be used for downstream analyses including:

UMAP visualization of cellular heterogeneity across WT and FAD samples.

Quantification of cell populations to assess changes in cell-type composition.

Heatmaps and UMAP overlays of marker genes to identify and spatially localize key cell types involved in Alzheimer’s pathology.

This pipeline lays the foundation for a comparative analysis of cell type-specific transcriptomic changes between healthy and Alzheimer’s-affected brain tissue.
