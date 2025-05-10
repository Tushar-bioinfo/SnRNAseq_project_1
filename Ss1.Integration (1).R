library(conflicted)
library(Seurat)
library(patchwork)
library(future)

options(Seurat.object.assay.version = "v3")
plan("multisession",workers=8)
options(future.globals.maxSize = 8 * 1024^3)
#Load 10x Matrices

s1.data <- Read10X(data.dir = "/Users/tusharsingh/Downloads/GSE227157_RAW/WT")
s2.data <- Read10X(data.dir = "/Users/tusharsingh/Downloads/GSE227157_RAW/FAD")

# Initialize the Seurat object with the raw data strained.
s1 <- CreateSeuratObject(counts = s1.data, project = "wt", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = s2.data, project = "fad", min.cells = 3, min.features = 200)

#mitochondrial percent is mt and ribosomal content is rp
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
#s1[["percent.rb"]] <- PercentageFeatureSet(s1, pattern = "^Rp[sl]")
s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^mt-")
#s2[["percent.rb"]] <- PercentageFeatureSet(s2, pattern = "^Rp[sl]")

#normalization of the data
#s1 <- NormalizeData(s1)
#s2<-NormalizeData(s2)

#scale data
all.genes.s1 <- rownames(s1)
s1_scale <-ScaleData(s1, features = all.genes.s1)
all.genes.s2 <- rownames(s2)
s2_scale <-ScaleData(s2, features = all.genes.s2)

#add metadata
s1_scale@meta.data[,"group"]<- "S1"
s2_scale@meta.data[,"group"]<- "S2"
gc()
merged_object <- merge(s1_scale, y = c(s2_scale), add.cell.ids = c("S1","S2"), project = "group")
treatment.list <- SplitObject(merged_object, split.by = "group")
#normalization and selection of  2000 variable features 

ifnb.list <- lapply(X = treatment.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize") 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

features <- SelectIntegrationFeatures(object.list = ifnb.list)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE) 
  x <- RunPCA(x, features = features, verbose = FALSE)})

muscle.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
muscle.combined <- IntegrateData(anchorset = muscle.anchors)

DefaultAssay(muscle.combined) <- "integrated"

#top10 <- head(VariableFeatures(muscle.combined), 10)

rm(s1.data)
