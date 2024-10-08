library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
sc_glucocort.data <- Read10X(data.dir = "../SCAT_CTL/")
# Initialize the Seurat object with the raw (non-normalized data).
sc_glucocort <- CreateSeuratObject(counts = sc_glucocort.data, project = "sc_glucocort", min.cells = 3, min.features = 200)
sc_glucocort
# Lets examine a few genes in the first thirty cells
sc_glucocort.data[c("Icam1", "Dpp4", "Limch1"), 1:30]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sc_glucocort[["percent.mt"]] <- PercentageFeatureSet(sc_glucocort, pattern = "^mt-")
# Show QC metrics for the first 5 cells
head(sc_glucocort@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(sc_glucocort, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(sc_glucocort, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_glucocort, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
sc_glucocort <- subset(sc_glucocort, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
#Normalization
sc_glucocort <- NormalizeData(sc_glucocort, normalization.method = "LogNormalize", scale.factor = 10000)
sc_glucocort <- NormalizeData(sc_glucocort)
# Identify variable features
sc_glucocort <- FindVariableFeatures(sc_glucocort, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_glucocort), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc_glucocort)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# Scaling the data
all.genes <- rownames(sc_glucocort)
sc_glucocort <- ScaleData(sc_glucocort, features = all.genes)
# example to remove unwanted variation
sc_glucocort <- ScaleData(sc_glucocort, vars.to.regress = "percent.mt")
# linear dimensional reduction
sc_glucocort <- RunPCA(sc_glucocort, features = VariableFeatures(object = sc_glucocort))
# Examine and visualize PCA results a few different ways
print(sc_glucocort[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc_glucocort, dims = 1:2, reduction = "pca")
DimPlot(sc_glucocort, reduction = "pca") + NoLegend()
DimHeatmap(sc_glucocort, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc_glucocort, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(sc_glucocort)
sc_glucocort <- FindNeighbors(sc_glucocort, dims = 1:19)
sc_glucocort <- FindClusters(sc_glucocort, resolution = 0.3)
# Look at cluster IDs of the first 5 cells
head(Idents(sc_glucocort), 5)

sc_glucocort <- RunUMAP(sc_glucocort, dims = 1:19)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(sc_glucocort, reduction = "umap", label = "TRUE")

saveRDS(sc_glucocort, file = "sc_glucocort.rds")
# find all markers of cluster 2
cluster2.markers <- FindMarkers(sc_glucocort, ident.1 = 2)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(sc_glucocort, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sc_glucocort.markers <- FindAllMarkers(sc_glucocort, only.pos = TRUE)
sc_glucocort.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

VlnPlot(sc_glucocort, features = c("Limch1", "F3", "Icam1", "Cd14", "Cd11b"))
# you can plot raw counts as well
VlnPlot(sc_glucocort, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(sc_glucocort, features = c("Mmp3", "Icam1", "Pdgfra"))

sc_glucocort.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(sc_glucocort, features = top10$gene) + NoLegend()
#Assigning cell type identity to clusters
new.cluster.ids <- c("EndoT1", "T-Cells", "EndoT2", "EndoT3", "ASC-ICAM1+", "ASC-Pi16+", "PeriC/Smooth",
                     "Schwann", "NK", "Macro", "Neutro", "B-Cells", "Neuron", "Erythroid")
names(new.cluster.ids) <- levels(sc_glucocort)
sc_glucocort <- RenameIdents(sc_glucocort, new.cluster.ids)
DimPlot(sc_glucocort, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
