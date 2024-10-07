library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
sc_glucocort.data <- Read10X(data.dir = "../SCAT_CTL/")
# Initialize the Seurat object with the raw (non-normalized data).
sc_glucocort <- CreateSeuratObject(counts = sc_glucocort.data, project = "sc_glucocort", min.cells = 3, min.features = 200)
sc_glucocort
# Lets examine a few genes in the first thirty cells
sc_glucocort.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
