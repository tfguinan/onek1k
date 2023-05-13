q()
# Not yet to be run
#################### This block stolen from _data for temporary testing ##########

library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(sctransform)
library(ggplot2)
library(pandoc)
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


features.average <- Connect("/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_features-average.h5seurat")
features.average <- LoadH5Seurat(features.average)
q()

########## All above is temporary ##########

library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(plotly)

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

## 1xall volcano plot
# LogFC by P
# load markers from RDS

# Load markers alongside log gene sums

logsums <- readRDS("/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_gene_log-sums.RDS")
old.dir <- getwd()
setwd('/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/1xall/')
files.list <- list.files(pattern = '\\.RDS$')

# Must add logsums at this stage
markers.list <- list()
for (file in files.list){
    file.markers <- readRDS(file)
    file.markers$gene <- rownames(file.markers)
    file.markers$log_sum <- logsums[file.markers$gene] 
    markers.list[[gsub('.RDS','',file)]] <- file.markers
}

# Markers will now have logsum and gene columns (not just rownames)
markers <- do.call(rbind, markers.list)
# rm(file.markers)
# rm(markers.list)
setwd(old.dir)


# average.matrix <- GetAssayData(features.average, slot = 'scale.data')
average.matrix <- GetAssayData(features.average, slot = 'counts')
colnames(average.matrix) <- gsub(' ', '_', colnames(average.matrix))
markers.list <- list()
old.dir <- getwd()
setwd('/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/1xall/')
files.list <- list.files(pattern = '\\.RDS$')
for (file in files.list){
    type <- strsplit(file, '-')[[1]][1]
    message(type)
    file.markers <- readRDS(file)
    file.markers$gene <- rownames(file.markers)
    file.markers$avg <- average.matrix[file.markers$gene, type]
    markers.list[[gsub('.RDS','', file)]] <- file.markers
} 
markers <- do.call(rbind, markers.list)


volcano_1 <- ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val), col=log_sum)) +
    geom_point()

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_volcano1_log.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
volcano_1
dev.off()

volcano_2 <- ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=log_sum)) +
    geom_point()

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_volcano2_log.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
volcano_2
dev.off()



## Do all 1x1 markers
# CD4 vs CD4

x1.markers <- readRDS('/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/1x1/CD4_Naive-v-CD4_CTL-1x1-lr.RDS')

volcano_3 <- ggplot(data=x1.markers, aes(x=avg_log2FC, y=-log10(p_val))) +
    geom_point()

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_volcano3_log.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
volcano_3
dev.off()

volcano_4 <- ggplot(data=x1.markers, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
    geom_point()

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_volcano4_log.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
volcano_3
dev.off()


## MA plot
# log2FC x average_expression
# See https://github.com/satijalab/seurat/issues/4167

features.average




ma_plot_1 <- ggplot(data = markers, aes(x=log_sum, y=avg_log2FC, col=p_val_adj)) + 
    geom_point(position = 'jitter')

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_ma_plot_1_jitter.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
ma_plot_1
dev.off()


ma_plot_2 <- ggplot(data = markers, aes(x=avg, y=avg_log2FC, col=p_val_adj)) + 
    geom_point(position = 'jitter')

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_ma_plot_2.pdf', width=34/2.54, height=20/2.54)                                                                                                                 
ma_plot_2
dev.off()

ma_plot_3 <- ggplot(data = markers, aes(x=avg, y=avg_log2FC, col=p_val_adj < 0.5)) + 
    geom_point() #+
    # scale_color_gradient(low = 'gray', high = 'red', breaks = c(FALSE, TRUE),
    #                    labels = c('p_val_adj >= 0.1', 'p_val_adj < 0.1'))

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/test_ma_plot_2.pdf', width=68/2.54, height=40/2.54)                                                                                                                 
ma_plot_3
dev.off()

################

features.average <- Connect("/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_features-average.h5seurat")
features.average <- LoadH5Seurat(features.average)

average.matrix <- GetAssayData(features.average, slot = 'counts')
scale.average.matrix <- GetAssayData(features.average, slot = 'scale.data')
log.average.matrix <- GetAssayData(features.average, slot = 'data')


matrix.list <- list(average.matrix, log.average.matrix, scale.average.matrix)
names(matrix.list) <- c('average', 'log.average', 'scale.average')
matrix.list <- lapply(matrix.list, function(x){colnames(x) <- gsub(' ', '_', colnames(x)); return (x)})

markers.list <- list()
old.dir <- getwd()
setwd('/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/1xall/')
files.list <- list.files(pattern = '\\.RDS$')
for (file in files.list){
    type <- strsplit(file, '-')[[1]][1]
    message(type)
    file.markers <- readRDS(file)
    file.markers$gene <- rownames(file.markers)
    for (name in names(matrix.list)){
        file.markers[name] <- matrix.list[[name]][file.markers$gene, type]
    }
    file.markers$mean.average <- mean.average.matrix
    markers.list[[gsub('.RDS','', file)]] <- file.markers
} 
rm(name)
markers <- do.call(rbind, markers.list)

mean.average.matrix <- rowMeans(average.matrix)

# markers$gene <- sapply(strsplit(rownames(markers), '\\.'), '[', 2)
markers$type <- sapply(strsplit(rownames(markers), '\\-'), '[', 1)

markers.p_adj.filtered <- subset(markers, p_val < 0.5)

# TODO add mean average (by gene)


# pct.1 - Percentage of cells where gene is detected in first group
# pct.2 - Percentage of cells where gene is detected in all others

plotly.p_adj.filtered <- plot_ly(data=markers.p_adj.filtered,
    x=~average, y=~avg_log2FC, 
    text=~paste('</br> gene: ', gene,
        '</br> type: ', type,
        '</br> p_val: ', p_val,
        '</br> p_val_adj: ', p_val_adj,
        '</br> %gene_cell_type: ', pct.1,
        '</br> %gene_all: ', pct.2))

# Unable to have selfcontained TRUE due to pandoc issues
htmlwidgets::saveWidget(
    widget = plotly.p_adj.filtered,
    file = "/data/menzies_projects/onek1k/share/TG/git_analysis/ma-p_avg-filtered.html",
    selfcontained=FALSE)


plotly.log.average <- plot_ly(data=markers, x=~log.average, y=~avg_log2FC, text=~name)

# Unable to have selfcontained TRUE due to pandoc issues
htmlwidgets::saveWidget(
    widget = plotly.log.average,
    file = "/data/menzies_projects/onek1k/share/TG/git_analysis/ma_log_average.html",
    selfcontained=FALSE)


plotly.scale.average <- plot_ly(data=markers, x=~scale.average, y=~avg_log2FC, text=~name)

# Unable to have selfcontained TRUE due to pandoc issues
htmlwidgets::saveWidget(
    widget = plotly.scale.average,
    file = "/data/menzies_projects/onek1k/share/TG/git_analysis/ma_scale_average.html",
    selfcontained=FALSE)

ma_plot <- ggplot(data = markers, aes(x=average, y=avg_log2FC, col=p_val_adj)) + 
    geom_point() +
    theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))


log_ma_plot <- ggplot(data = markers, aes(x=log.average, y=avg_log2FC, col=p_val_adj)) + 
    geom_point() +
    theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

scaled_ma_plot <- ggplot(data = markers, aes(x=scale.average, y=avg_log2FC, col=p_val_adj)) + 
    geom_point() +
    theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))


pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/ma_plot.pdf', width=68/2.54, height=40/2.54)     
ma_plot
dev.off()

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/ma_plot_log.pdf', width=68/2.54, height=40/2.54)
log_ma_plot
dev.off()
    
pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/ma_plot_scaled.pdf', width=68/2.54, height=40/2.54)
scaled_ma_plot
dev.off()


library(ComplexHeatmap)


cell_type_regulon <- read.csv('/data/menzies_projects/onek1k/share/TG/git_analysis/cell_type_rss.csv')

pdf('/data/menzies_projects/onek1k/share/TG/git_analysis/cell_type-onek1k-heatmap.pdf', width=41.5/2.54, height=68.5/2.54)                                                                                                                 
Heatmap(cell_type_regulon,
    row_names_gp = grid::gpar(fontsize = 2), column_names_gp = grid::gpar(fontsize = 15))
dev.off()


# Load necessary libraries
library(ggplot2)
library(reshape2)

# Read in CSV file
df <- read.csv('/data/menzies_projects/onek1k/share/TG/git_analysis/cell_type_rss.csv', header=TRUE, row.names=1)




library(tidyverse)

dt2 <- df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)


# Create heatmap
ggplot(dt2, aes(x=rowname, y=colname, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient(low="white", high="red")

# Save plot as PDF
ggsave('/data/menzies_projects/onek1k/share/TG/git_analysis/cell_type-onek1k-heatmap.pdf')
