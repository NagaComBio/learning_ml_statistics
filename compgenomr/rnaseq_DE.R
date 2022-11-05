## DE analysis
library(DESeq2)
library(ggplot2)
counts_file <-  system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                            package='compGenomRData')
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package='compGenomRData')

counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
counts <- as.matrix(subset(counts, select = c(-width)))
colData <- read.table(coldata_file, header = T, sep = '\t', stringsAsFactors = T)

d_formule <- '~ group'

## Create deseq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = as.formula(d_formule))

dds

## remove genes with low expression
dds <- dds[rowSums(counts(dds)) > 1,]

## perform DE
dds <- DESeq(dds)

DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
DEresults <- DEresults[order(DEresults$pvalue),]
print(DEresults)
## diagnostic plot
plotMA(dds)

ggplot(data = as.data.frame(DEresults),
       aes(x = pvalue)) + 
  geom_histogram(bins = 100)


library(EnhancedVolcano)
EnhancedVolcano(DEresults, x='log2FoldChange', y = 'padj', lab = rownames(DEresults))

######## Functional enrichment
library(gProfileR)
library(knitr)

DE <- DEresults[!is.na(DEresults$padj),]
DE <- DE[DE$padj < 0.1,]
DE <- DE[abs(DE$log2FoldChange) > 1,]
dim(DEresults)
dim(DE)

goi <- rownames(DE)
goResults <- gprofiler(query = goi,
                       organism = 'hsapiens',
                       src_filter = 'GO',
                       hier_filtering = 'moderate')
head(goResults)

## GSEA
goResults <- goResults[order(goResults$p.value),]
go <- goResults[goResults$overlap.size < 100,]
geneset1 <- unlist(strsplit(go[1,]$intersection, ','))

normalzied_counts <- DESeq2::counts(dds, normalized = T)
geneset2 <- sample(rownames(normalzied_counts), 25)

geneSets <- list('top_GO_term' = geneset1,
                 'random_set' = geneset2)

library(gage)
gseaResults <- gage(exprs = log2(normalzied_counts+1),
                    ref = match(rownames(colData[colData$group == 'CTRL',]),
                                colnames(normalzied_counts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalzied_counts)),
                    gsets = geneSets, compare = 'as.group')

gseaResults

## heatmap
library(pheatmap)
M <- normalzied_counts[rownames(normalzied_counts) %in% geneset1, ]
pheatmap(log2(M+1),
         annotation_col = colData,
         show_rownames = T,
         fontsize = 8,
         scale = 'row',
         cutree_rows = 2,
         cutree_cols = 2)
