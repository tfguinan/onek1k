library(Seurat)
library(SeuratDisk)
library(ggplot2)


# Wrapper for Seurat VlnPlot, allowing plotting of individuals and pools
v.plot <- function(object, genes, ids, pools, assay, slot, path, name){
    Idents(object) <- 'individual'
    for (gene in genes){
        plot <- VlnPlot(object=object, features=gene, idents=ids, assay=assay, slot=slot) + 
            theme(legend.position='none')
        pdf(sprintf('%s/%s-%s-individual-violin-%s-%s.pdf', path, name, gene, assay, slot))
        print(plot)
        dev.off()
    }

    Idents(object) <- 'pool'
    for (gene in genes){
        plot <- VlnPlot(object=object, features=gene, idents=pools, assay=assay, slot=slot) + 
            theme(legend.position='none')
        pdf(sprintf('%s/%s-%s-pool-violin-%s-%s.pdf', path, name, gene, assay, slot))
        print(plot)
        dev.off()    
    }
}


h5.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/seurat/onek1k_main.h5seurat'
h5.file <- SeuratDisk::Connect(h5.path)
onek1k.main <- SeuratDisk::LoadH5Seurat(h5.file, assays = c("RNA", "SCT"), reductions = NULL, graphs = FALSE, images = FALSE)
h5.file$close_all()

message('Dataset head:')
head(onek1k.main)

Idents(onek1k.main) <- 'predicted.celltype.l2'

# TODO improve subsetting to enable a list of objects
onek1k.plasmablast <- subset(x = onek1k.main, idents = 'Plasmablast')
onek1k.mono <- subset(x = onek1k.main, idents = c('CD14 Mono', 'CD16 Mono'))

Idents(onek1k.plasmablast) <- 'individual'
Idents(onek1k.mono) <- 'individual'

message('Subset head:')
head(onek1k.plasmablast)
message('Subset head:')
head(onek1k.mono)

plot.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis'
# plot.path <- '/data/menzies_projects/onek1k/share/TG/git_analysis/test_visualisation'

gene.c <- c('HSP90B1', 'MCL1', 'RPS27', 'TNFRSF17')
pool.c <- levels(onek1k.main$pool)


# Plot genes on UMAP
for (gene in gene.c){
    plot <- FeaturePlot(onek1k.main, features = gene)
    pdf(sprintf('%s/%s_feature.pdf', plot.path, gene), width=58.5/2.54, height=58.5/2.54)                                                                                                                 
    print(plot)
    dev.off()
}


data.plasmablast <- FetchData(onek1k.plasmablast, append(gene.c, 'ident'))
# Find individual IDs based on # cells
indv.count.plasmablast <- sort(table(data.plasmablast$ident), decreasing=TRUE)
indv.id.plasmablast <- names(indv.count.plasmablast)[1:20]

print(indv.id.plasmablast)

rm(data.plasmablast, indv.count.plasmablast)


# Plotting with wrapper function
v.plot(onek1k.plasmablast, gene.c, indv.id.plasmablast, pool.c, 'RNA', 'counts', plot.path, 'plasmablast')

v.plot(onek1k.plasmablast, gene.c, indv.id.plasmablast, pool.c, 'SCT', 'counts', plot.path, 'plasmablast')


data.mono <- FetchData(onek1k.mono, append(gene.c, 'ident'))
# Find individual IDs based on # cells
indv.count.mono <- sort(table(data.mono$ident), decreasing=TRUE)
indv.id.mono <- names(indv.count.mono)[1:20]

print(indv.id.mono)

rm(data.mono, indv.count.mono)


# Plotting with wrapper function
v.plot(onek1k.mono, gene.c, indv.id.mono, pool.c, 'RNA', 'counts', plot.path, 'mono')

v.plot(onek1k.mono, gene.c, indv.id.mono, pool.c, 'SCT', 'counts', plot.path, 'mono')


print('Done')
q()



# Old code for selecting individual IDs
# ### Individual mean; gene mean selection method
# # Aggregate (take mean) of all individuals cells
# data <- aggregate(. ~ ident, data = data, mean)

# # Sort by mean of all gene columns
# data$mean <- apply(data[,gene.c], 1, mean)
# data <- data[order(-data$mean),]

# rownames(data) <- NULL

# message('Individual ranking data head:')
# print(head(data))

# incr <- nrow(data)/10

# indv <- c(1,
#     1 + round(incr),
#     1 + round(2 * incr),
#     1 + round(3 * incr),
#     1 + round(4 * incr),
#     1 + round(5 * incr),
#     1 + round(6 * incr),
#     1 + round(7 * incr),
#     1 + round(8 * incr),
#     nrow(data)
# )

# indv.id <- as.character(data$ident[indv])
# rm(data, incr, indv)
# ###

# print(indv.id)

# # TODO look into DotPlots by celltype, individual
# Idents(onek1k.main) <- 'predicted.celltype.l2'

# plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c)
# pdf(sprintf('%s/%s_celltype_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
# plot
# dev.off()

# Idents(onek1k.main) <- 'individual'

# plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c, idents = indv.id)
# pdf(sprintf('%s/%s_individual_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
# plot
# dev.off()

# Idents(onek1k.main) <- 'pool'

# plot <- DotPlot(onek1k.main, assay = 'SCT', features = gene.c)
# pdf(sprintf('%s/%s_pool_dotplot_sct_data.pdf', plot.path, gene), width=34/2.54, height=20/2.54)                                                                                                                 
# plot
# dev.off()

