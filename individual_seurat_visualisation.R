# Plot violin for 10 individuals on a single gene
library(Seurat)
library(SeuratDisk)
library(ggplot2)

h5.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_main.h5seurat'
h5.file <- SeuratDisk::Connect(h5.path)
onek1k.main <- SeuratDisk::LoadH5Seurat(h5.file, assays = c("RNA", "SCT"), reductions = NULL, graphs = FALSE, images = FALSE)
h5.file$close_all()

plot.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis'

# Plot normalized expression
DefaultAssay(onek1k.main) <- 'SCT'
Idents(onek1k.main) <- 'individual'


# This block finds individual IDs between min and max
gene.c <- c('HSP90B1', 'MCL1', 'RPS27')

data <- FetchData(onek1k.main, append(gene.c, 'ident'))

data <- aggregate(. ~ ident, data = data, mean)

# Sort by mean of all gene columns
data$mean <- apply(data[,gene.c], 1, mean)
data <- data[order(-data$mean),]

rownames(data) <- NULL

incr <- nrow(data)/10

indv <- c(1,
    1 + round(incr),
    1 + (2 * round(incr)),
    1 + (3 * round(incr)),
    1 + (4 * round(incr)),
    1 + (5 * round(incr)),
    1 + (6 * round(incr)),
    1 + (7 * round(incr)),
    1 + (8 * round(incr)),
    nrow(data)
)

indv.id <- as.character(data$ident[indv])

rm(data, incr, data, gene.c)


for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_sct.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    plot
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_sct.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    plot
    dev.off()
}


# Plot genes in UMAP
for (gene in gene.c){
    plot <- FeaturePlot(onek1k.main, features = gene)
    pdf(sprintf('%s/%s_feature.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    plot
    dev.off()
}


# Plot unnormalized expression
DefaultAssay(onek1k.main) <- 'RNA'
Idents(onek1k.main) <- 'individual'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene, idents = indv.id) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_individual_violin_sct.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    plot
    dev.off()
}

Idents(onek1k.main) <- 'pool'

for (gene in gene.c){
    plot <- VlnPlot(object = onek1k.main, features = gene) + theme(legend.position = 'none')
    pdf(sprintf('%s/%s_pool_violin_sct.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
    plot
    dev.off()
}