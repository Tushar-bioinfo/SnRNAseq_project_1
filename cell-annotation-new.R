library(Seurat)
library(patchwork)
library(dplyr)
library(clustree)
# Process integrated dataset, 10 PCs
treatment.integrated_scale <- ScaleData(muscle.combined)
treatment.integrated_scale <- RunPCA(treatment.integrated_scale, 
                                     npcs = 10,
                                     approx = FALSE)
treatment.integrated_scale <- RunUMAP(treatment.integrated_scale,reduction = "pca", features = NULL, 
                                      dims = 1:10)
DimPlot(treatment.integrated_scale, reduction = "umap")

treatment.integrated_scale <- FindNeighbors(treatment.integrated_scale, 
                                            dims = 1:10)


#Determine Clustering resolution
treatment.integrated_scale_cluster <- FindClusters(treatment.integrated_scale, 
                                                   resolution = c(0,0.1,0.2,0.3,0.4,0.5))
View(treatment.integrated_scale_cluster@meta.data)

clustree(treatment.integrated_scale_cluster, prefix = "integrated_snn_res.")
#Choose Clustering resolution
treatment.integrated_scale <-FindClusters(treatment.integrated_scale, resolution = 0.3)
DimPlot(treatment.integrated_scale, reduction = "umap",label=T)
#rm(treatment.integrated_scale_cluster)
#Find conserved markers for every cluster
conserved.marker.0 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="0")
conserved.marker.1 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="1")
conserved.marker.2 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="2")
conserved.marker.3 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="3")
conserved.marker.4 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="4")
conserved.marker.5 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group",logfc.threshold=1,
                                           ident.1="5")
conserved.marker.6 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="6")
conserved.marker.7 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="7")
conserved.marker.8 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="8")
conserved.marker.9 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                           ident.1="9")
conserved.marker.10 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="10")
conserved.marker.11 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="11")
conserved.marker.12 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="12")
conserved.marker.13 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="13")
conserved.marker.14 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="14")
conserved.marker.15 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="15")
conserved.marker.16 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="16")
conserved.marker.17 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
                                            ident.1="17")
#conserved.marker.18 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
 #                                           ident.1="18")
#conserved.marker.19 <- FindConservedMarkers(treatment.integrated_scale, only.pos=TRUE, grouping.var = "group", logfc.threshold=1,
 #                                           ident.1="19")
#Form new Seurat object with clusters expressing strong conserved markers
ywt_fad <- subset(treatment.integrated_scale, idents = c(0,1,2,3,5,6,7,8,10,11,13,14,15,16,17))
#Annotate Cell-types in order of subset idents
new.cluster.ids <- c("Microglia","Oligodendrocyte","Astrocyte","Endothelial","Oligodendrocyte","Neuron",
                     "OPC","Oligodendrocyte","Oligodendrocyte","OPC","Oligodendrocyte","Neuron","Endothelial",
                     "Oligodendrocyte","Oligodendrocyte")
names(new.cluster.ids) <- levels(ywt_fad)
ywt_fad <- RenameIdents(ywt_fad, new.cluster.ids)

#Place Cell-types in alphabetical order
levels(ywt_fad) <- c("Astrocyte","Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC")

#cell type annotation
DimPlot(ywt_fad, reduction = "umap")

DimPlot(treatment.integrated_scale, reduction = "umap", group.by = "group", label = TRUE) + ggtitle("UMAP: WT vs FAD (Unannotated)")
DimPlot(ywt_fad, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP: Annotated Cell Types")



# Extract metadata and cell identities
meta_data <- ywt_fad@meta.data
meta_data$celltype <- Idents(ywt_fad)
meta_data$group <- recode(meta_data$group, "S1" = "WT", "S2" = "FAD")
# Count cells per group and cell type
cell_counts_df <- meta_data %>%
  group_by(group, celltype) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count)) %>%
  ungroup()
print(cell_counts_df)

# Plot the proportions of each cluster
library(ggplot2)

ggplot(cell_counts_df, aes(x = group, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1)),
            position = position_fill(vjust = 0.5), size = 3, color = "black") +
  ylab("Proportion of Cells") +
  xlab("Group") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  ggtitle("Cell Type Proportions in WT and FAD Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Heatmap for marker genes
marker.genes <- c('Slc1a2', 'Plpp3','Clu','Flt1', 'Cldn5','Ly6a','Btg2', 'Kctd12', 
                  'Cx3cr1', 'Gria2', 'Nrxn1', 'Syne1', 'Enpp2', 'Ptgds' ,'Fabp7', 'Apod', 'Ptn')
DefaultAssay(ywt_fad) <- "RNA"

ywt_fad <- ScaleData(ywt_fad, features = marker.genes)

DoHeatmap(ywt_fad, features = marker.genes, group.by = "ident") +
  ggtitle("Heatmap of Marker Gene Expression")

FeaturePlot(ywt_fad, features = marker.genes, 
            cols = c("lightgrey", "blue"), reduction = "umap", 
            pt.size = 0.4, ncol = 6)

#VlnPlot(ywt_fad, features = c("Gfap", "Aif1", "Pdgfra", "Snap25"), 
       # group.by = "orig.ident", pt.size = 0)
