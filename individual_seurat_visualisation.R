# Plot violin for 10 individuals on a single gene
# Plot violin for all pools on single gene
# Plot dotplot for 10 individuals on gene list
# Plot dotplot for all pools on gene list
library(Seurat)
library(SeuratDisk)
library(ggplot2)


h5.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_main.h5seurat'
# h5.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/celltype_split/NK-onek1k.h5seurat'
h5.file <- SeuratDisk::Connect(h5.path)
onek1k.main <- SeuratDisk::LoadH5Seurat(h5.file, assays = c("RNA", "SCT"), reductions = NULL, graphs = FALSE, images = FALSE)
h5.file$close_all()

message('Dataset head:')
head(onek1k.main)

plot.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis'
# plot.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/test_visualisation'

DefaultAssay(onek1k.main) <- 'SCT'
Idents(onek1k.main) <- 'individual'


# This block finds individual IDs between min and max
gene.c <- c('HSP90B1', 'MCL1', 'RPS27', 'TNFRSF17')
# Need to match syntax for selecting pools

# pool.c <- c('pool_1', 'pool_2', 'pool_3', 'pool_4', 'pool_5',
#     'pool_6', 'pool_7', 'pool_8', 'pool_9', 'pool_10')

pool.c <- levels(onek1k.main$pool)


data <- FetchData(onek1k.main, append(gene.c, 'ident'))

data <- aggregate(. ~ ident, data = data, mean)

# Sort by mean of all gene columns
data$mean <- apply(data[,gene.c], 1, mean)
data <- data[order(-data$mean),]

rownames(data) <- NULL

message('Individual ranking data head:')
print(head(data))

incr <- nrow(data)/10

indv <- c(1,
    1 + round(incr),
    1 + round(2 * incr),
    1 + round(3 * incr),
    1 + round(4 * incr),
    1 + round(5 * incr),
    1 + round(6 * incr),
    1 + round(7 * incr),
    1 + round(8 * incr),
    nrow(data)
)

indv.id <- as.character(data$ident[indv])
print(indv.id)

rm(data, incr, indv)

# TODO look at KDE adjustment

# Plot normalized expression
DefaultAssay(onek1k.main) <- 'SCT'

### Plot SCT corrected counts
Idents(onek1k.main) <- 'individual'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id,
        assay = "SCT", slot = "counts") + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_sct_counts.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = pool.c,
        assay = "SCT", slot = "counts", pt.size = 0) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_sct_counts.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

### Plot SCT corrected log-normalised counts
Idents(onek1k.main) <- 'individual'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id,
        assay = "SCT", slot = "data") + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = pool.c,
        assay = "SCT", slot = "data", pt.size = 0) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}


### Plot genes on UMAP
for (gene in gene.c){
    plot <- FeaturePlot(onek1k.main, features = gene)
    pdf(sprintf('%s/%s_feature.pdf', plot.path, gene), width=58.5/2.54, height=58.5/2.54)                                                                                                                 
    print(plot)
    dev.off()
}


# Plot unnormalised expression
DefaultAssay(onek1k.main) <- 'RNA'

### Plot RNA raw counts
Idents(onek1k.main) <- 'individual'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id,
        assay = "RNA", slot = "counts") + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_rna_counts.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = pool.c,
        assay = "RNA", slot = "counts", pt.size = 0) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_rna_counts.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

### Plot RNA log-normalised counts
Idents(onek1k.main) <- 'individual'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id,
        assay = "RNA", slot = "data") + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_rna_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = pool.c, 
        assay = "RNA", slot = "data", pt.size = 0) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_rna_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    print(plot)
    dev.off()
}


# TODO look into DotPlots by celltype, individual


Idents(onek1k.main) <- 'predicted.celltype.l2'

plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c)
pdf(sprintf('%s/%s_celltype_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
plot
dev.off()

Idents(onek1k.main) <- 'individual'

plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c, idents = indv.id)
pdf(sprintf('%s/%s_individual_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
plot
dev.off()

Idents(onek1k.main) <- 'pool'

plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c)
pdf(sprintf('%s/%s_pool_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
plot
dev.off()