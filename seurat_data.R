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
rm.list <- c('Platelet', 'Eryth', 'Doublet')
onek1k_main <- SeuratObject:::subset.Seurat(onek1k_main, idents=rm.list, invert=TRUE) 
print(levels(onek1k_main))


# Use celltype pairing metadata to run FindMarkers 
# pairs.df <- read.csv(path(raw.path, 'full_pairs.csv'))
pairs.df <- read.csv(pairs.path)

for (x in 1:nrow(pairs.df)){
    celltype.1 <- pairs.df$celltype.1[x]
    celltype.2 <- pairs.df$celltype.2[x]

    # Using predicted.celltype.l2:
    # Create file path friendly string for celltype pairing
    celltype.pair <- paste(gsub(' ', '_', celltype.1), gsub(' ', '_', celltype.2), sep = celltype.sep)
    message(celltype.pair)

    rds.path <- path(v1.path, sprintf("%s-%s-%s.RDS", gsub(' ', '_', celltype.pair), '1x1', 'lr'))
    # e.g. /project_directory_base_path/output_data/seurat/1v1/CD4_TEM-x-CD4_TCM-1x1-lr.RDS
    message(rds.path)
    if (file.exists(rds.path)){
        message(sprintf('This path is preexisting: %s', rds.path))
        # markers <- readRDS(rds.path) 
    } else {
        message(sprintf('Running FindMarkers for %s', celltype.pair))
        
        tryCatch(
        {markers <- Seurat::FindMarkers(object = onek1k_main,
            ident.1 = celltype.1,
            ident.2 = celltype.2,
            min.pct = 0.1,
            test.use = 'LR')
            message(sprintf('Writing to: %s', rds.path))
            saveRDS(markers, rds.path)
            }, error = function(e) {print('warning: likely missing celltype (or wrong ident)')})

    }
}
# TODO rerun markers with a different test?

# Save object with feature average by celltype
SeuratObject::DefaultAssay(onek1k_main) <- "RNA"
# Equal matrix in both slots indicates no log normalisation done
# identical(GetAssayData(onek1k_main, assay = 'RNA', slot = 'data'), GetAssayData(onek1k_main, assay = 'RNA', slot = 'counts'))

# Use counts slot as input for AverageExpression
# Will perform ScaleData by default
# See https://satijalab.org/seurat/reference/averageexpression 
features.average <- Seurat::AverageExpression(object = onek1k_main, slot = 'counts', return.seurat = TRUE)
SeuratDisk::SaveH5Seurat(features.average, filename = path(out.path, sprintf('%s_features-average.h5seurat', data.name)))

DefaultAssay(onek1k_main) <- "SCT"
for (focus.celltype in levels(onek1k_main)){
    message(focus.celltype)

    rds.path <- path(vall.path, sprintf("%s-%s-%s.RDS", gsub(' ', '_', focus.celltype), '1xall', 'lr'))
    # e.g. /project_directory/output_data/seurat/1vall/CD4_TCM-1xall-lr.RDS
    message(rds.path)
    if (file.exists(rds.path)){
        message(sprintf('This path is preexisting: %s', rds.path))
        # markers <- readRDS(rds.path) 
    } else {
        message(sprintf('Running FindMarkers for %s', focus.celltype))
        
        tryCatch(
        {markers <- Seurat::FindMarkers(object = onek1k_main,
            ident.1 = focus.celltype,
            min.pct = 0.1,
            test.use = 'LR')
            message(sprintf('Writing to: %s', rds.path))
            saveRDS(markers, rds.path)
            }, error = function(e) {print('warning: likely missing celltype (or wrong ident)')})
    }
}

# TODO implement celltype_split here:
# Celltype split with SCT assay 
SeuratDisk::SaveH5Seurat(onek1k_main, filename = path(out.path, sprintf('%s_main.h5seurat', data.name)))
saveRDS(onek1k_main, path(out.path, sprintf('%s_main.RDS', data.name)))
# Convert h5seurat to h5ad, for import to Python
# TODO look into in place conversion
# See https://mojaveazure.github.io/seurat-disk/reference/Convert.html
SeuratDisk::Convert(path(out.path, sprintf('%s_main.h5seurat', data.name)), dest = 'h5ad')
rm(onek1k_main)