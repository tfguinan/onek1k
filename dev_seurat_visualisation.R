########### This block stolen from _data for temporary testing ##########

library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(sctransform)
library(ggplot2)
# Seems dataset is too large for parallelising (duplicate object in memory?)
# library(future)


# TODO check all functions are used +1 time, else use code directly

path <- function(path.1, path.2, create = FALSE){
    path.out <- paste(path.1, path.2, sep = "/")
    if (create){
        if (!file.exists(path.out)){
            dir.create(path.out)
            # see mkdir -p
            message('Created ', path.out)
        } else {message(path.out, ' already exists')}
    }
    return(path.out)
}


# Returns a generated path for an RDS file containing cell, analysis, and test information
build.rds_path <- function(path.in, focus.celltype, analysis.type, test.used){
    return(path(path.in, sprintf("%s-%s-%s.RDS", gsub(' ', '_', focus.celltype), analysis.type, test.used)))
}


args <- commandArgs(trailingOnly=TRUE)
argslen <- length(args)

if (argslen < 4) stop('error: insufficient arguments (< 4)')

# # Temporary
# base.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis'
# data.name <- 'onek1k'
# h5.path <- '/data/menzies_projects/onek1k/share/share_raw/onek1k_seurat_v210819.h5seurat'
# pairs.path <- '/data/menzies_projects/onek1k/share/share_raw/full_pairs.csv'
# # Optionally add TOB ID to h5seurat object (in memory)
# id.path <- '/data/menzies_projects/onek1k/share/share_raw/OneK1K_sample_IDs_210304.tsv'

base.path <- args[1]
data.name <- args[2]
h5.path <- args[3]
pairs.path <- args[4]
# Optionally add TOB ID to h5seurat object (in memory)
id.path <- if (argslen < 5) NULL else args[5]

setwd(base.path)

raw.path <- path(base.path, 'raw_data', create = FALSE)
# Should fail (already exists from setup.sh)
out.path <- path(base.path, 'output_data/seurat', create = FALSE)


# The seperator used when joining cell types into a file name
celltype.sep <- '-v-'

# For RDS files of FindMarkers for 1x1/all
vall.path <- path(out.path, sprintf('1%sall', celltype.sep), create = TRUE)
v1.path <- path(out.path, sprintf('1%s1', celltype.sep), create = TRUE)
# For plots common to all analyses
# Part of initial setup (assumed)

h5.file <- SeuratDisk::Connect(h5.path)
onek1k_main <- SeuratDisk::LoadH5Seurat(h5.file, assays = c("RNA", "SCT"), reductions = NULL, graphs = FALSE, images = FALSE)

# Here we add additional required ID information
onek1k_main[['old.ident']] <- Seurat::Idents(object = onek1k_main)
SeuratObject::Idents(object = onek1k_main) <- 'individual'

# Full data is in memory, close file connection
h5.file$close_all()

if (!is.null(id.path)){
    id.csv <- read.csv(file = id.path, sep = '\t')

    # Create named vector
    id.map <- id.csv$TOB_Dark_ID
    names(id.map) <- id.csv$PERSON
    # Create empty data in Seurat object
    onek1k_main$TOB_Dark_ID <- NA
    # Map TOB to individual ID per cell
    onek1k_main$TOB_Dark_ID <- id.map[onek1k_main$individual]
    if (sum(is.na(onek1k_main$TOB_Dark_ID)) > 0){
        message("warning: TOB ID sum NA > 0")
        warning(paste(na.sum))
    }
}



# Return default identity (Idents) to cell type predictions
SeuratObject::Idents(object = onek1k_main) <- 'predicted.celltype.l2'
print(levels(onek1k_main))
# Subset main object to exclude Doublets, Erythrocytes and Platelets
# In place to conserve memory
rm.list <- c('Platelet', 'Eryth', 'Doublet')
onek1k_main <- SeuratObject:::subset.Seurat(onek1k_main, idents=rm.list, invert=TRUE) 
print(levels(onek1k_main))


# Use celltype pairing metadata to run FindMarkers 
# pairs.df <- read.csv(path(raw.path, 'full_pairs.csv'))
pairs.df <- read.csv(pairs.path)

message("Start is all ok")

DefaultAssay(onek1k_main) <- "RNA"
# Use counts slot from original h5seurat to calculate averages
features.average <- Seurat::AverageExpression(object = onek1k_main, slot = 'counts', return.seurat = TRUE)
SeuratDisk::SaveH5Seurat(features.average, filename = path(out.path, sprintf('%s_features-average.h5seurat', data.name)))


q()

########## All above is temporary ##########

library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)

# TODO
# TODO sanitise all cell type names, remove spaces

# Path handling
# Ident checks?
# Load full object (onek1k_main)

## QC (supp?)

pdf(path(plot.path, "ncount_RNA_vs_percent_MT.pdf"), width=58.5/2.54, height=58.5/2.54)
Seurat::FeatureScatter(onek1k_main, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

# PCA visualisation
pdf(path(plot.path, "VizDimLoadings1-pca.pdf"), width=58.5/2.54, height=58.5/2.54)
Seurat::VizDimLoadings(onek1k_main, dims = 1:25, reduction = "pca")
dev.off()

pdf(path(plot.path, "VizDimLoadings2-pca.pdf"), width=58.5/2.54, height=58.5/2.54)
Seurat::VizDimLoadings(onek1k_main, dims = 26:50, reduction = "pca")
dev.off()

### Figure 1 ###

## UMAP
pdf(path(plot.path, "Onek1k_UMAP.pdf"), width=58.5/2.54, height=58.5/2.54)
DimPlot(onek1k_main, reduction = 'umap',
 order = cell.list, cols = col.cells, shuffle = TRUE,
 label = TRUE, raster = TRUE, raster.dpi=c(2048, 2048)) + ggtitle('Onek1k UMAP')
dev.off()
# Also save projection coordinates (possibly done in seurat_data)

## Heatmap
# Uses complexheatmap package

# Need to source these markers from seurat_data
# Also no hardcoded paths!
upper.markers <- read.csv('/data/menzies_projects/onek1k/share/analysis/DGE/output/0-585-1xall_upper_markers.csv')
lower.markers <- read.csv('/data/menzies_projects/onek1k/share/analysis/DGE/output/0-585-1xall_lower_markers.csv')
markers <- rbind(upper.markers['gene'], lower.markers['gene'])

# Take character vector (x$y method)
unique.markers <- unique(markers)$gene

# Get average expression scaled data 
# ScaleData performed by default when counts slot is used as input (see features.average in seurat_data.R)
scale.data <- GetAssayData(features.average, assay = 'RNA', slot = 'scale.data')

# Default (complete clustering)
pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/features_average_heatmap_complete.pdf', width=41.5/2.54, height=68.5/2.54)                                                                                                                 
Heatmap(scale.data[unique.markers,],
 row_names_gp = grid::gpar(fontsize = 2), column_names_gp = grid::gpar(fontsize = 15))
dev.off()

# Improved column clustering (T cells)
pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/features_average_heatmap_r_complete_c_single.pdf', width=41.5/2.54, height=68.5/2.54)                                                                                                                 
Heatmap(scale.data[unique.markers,],
 row_names_gp = grid::gpar(fontsize = 2), column_names_gp = grid::gpar(fontsize = 15),
 clustering_method_rows = 'complete', clustering_method_columns = 'single')
dev.off()

# # Using ward clustering
# pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/features_average_heatmap_groupby_subset_ward_data.pdf', width=41.5/2.54, height=68.5/2.54)                                                                                                                 
# Heatmap(scale.data[unique.markers,],
#  row_names_gp = grid::gpar(fontsize = 2), column_names_gp = grid::gpar(fontsize = 15),
#  clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2')
# dev.off()

# Feature plots

### Figure 2 ###
## Matrix (TODO assess)

## 1v1 Violins