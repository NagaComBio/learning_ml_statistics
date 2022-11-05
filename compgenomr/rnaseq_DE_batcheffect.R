library(pheatmap)
library(DESeq2)

counts_file <- system.file('extdata/rna-seq/SRP021193.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP021193.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header=T, sep='\t',
                      stringsAsFactors = T)
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_vector')

genelength <- counts$width
rpk <- apply(subset(counts, select = c(-width)), 2, function(x){
  x / (genelength/1000)
})

tpm <- apply(rpk, 2,  FUN = function(x){
  (x / sum(x))/10**6
})

selected_genes <- names(sort(apply(tpm, 1, var), decreasing = T))[1:100]

pheatmap(tpm[selected_genes,],
         scale = 'row',annotation_col = colData)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(subset(counts, select = c(-width))),
                              colData = colData,
                              design = as.formula(~ LibrarySelection + group))

dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
DEresults <- DEresults[order(DEresults$pvalue),]

normalized_counts <- counts(dds, normalized = T)
de_genes <- rownames(DEresults)[1:500]
pheatmap(normalized_counts[de_genes,],
         scale = 'row', 
         annotation_col = colData, 
         show_rownames = F, cutree_cols = 5)

(is.na(apply(tpm[de_genes,], 1, scale)))


##########################
## Automatically estimate covariates
counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, 
                      sep = '\t', stringsAsFactors = TRUE)
# simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')

#find gene length normalized values 
geneLengths <- counts$width
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
selectedGenes <- names(sort(apply(tpm, 1, var), 
                            decreasing = T)[1:100])
pheatmap(tpm[selectedGenes,], 
         scale = 'row',
         annotation_col = colData, 
         cutree_cols = 2, 
         show_rownames = FALSE)

library(EDASeq)
# remove 'width' column from counts
countData <- as.matrix(subset(counts, select = c(-width)))
# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)

par(mfrow = c(1, 2))
plotRLE(set, outline=F, ylim = c(-4, 4), col = as.numeric(colData$group))
plotPCA(set, col=as.numeric(colData$group),
        adj = 0.5,
        ylim = c(-0.7, 0.5), 
        xlim = c(-0.5, 0.5))
### On tpm
plotRLE(tpm, outline=F, ylim = c(-4, 4), col = as.numeric(colData$group))
plotPCA(tpm, col=as.numeric(colData$group),
        adj = 0.5,
        ylim = c(-0.7, 1), 
        xlim = c(-0.5, 0.5))

## Reading in house keeping genes
library(RUVSeq)

hk_genes <- read.table(file = system.file("extdata/rna-seq/HK_genes.txt", 
                                          package = 'compGenomRData'), 
                       header = FALSE)

hk_genes <- intersect(rownames(set), hk_genes$V1)

par(mfrow = c(2,2))
for(k in 1:4){
  set_g <- RUVg(x = set, cIdx = hk_genes, k = k)
  plotPCA(set_g, col = as.numeric(colData$group), 
          cex = 0.9, adj = 0.5, 
          main = paste0('with RUVg, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1))
}

# choose k = 1
set_g <- RUVg(x = set, cIdx = hk_genes, k = 1)

# RLA plots
par(mfrow = c(1, 2))
plotRLE(set, outline = F, col=as.numeric(colData$group), ylim = c(-2, 2))
plotRLE(set_g, outline = F, col=as.numeric(colData$group), ylim = c(-2, 2))

# PCA plots
par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVg', 
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
plotPCA(set_g, col=as.numeric(colData$group), adj = 0.5, 
        main = 'with RUVg',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))

# RUVs
difference <- makeGroups(colData$group)
difference

par(mfrow = c(2, 2))
for(k in 1:4) {
  set_s <- RUVs(set, unique(rownames(set)), 
                k=k, difference) #all genes
  plotPCA(set_s, col=as.numeric(colData$group), 
          cex = 0.9, adj = 0.5, 
          main = paste0('with RUVs, k = ',k), 
          ylim = c(-1, 1), xlim = c(-0.6, 0.6))
}
set_s <- RUVs(set, unique(rownames(set)), k=3, difference) #
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), 
        main = 'without RUVs')
plotRLE(set_s, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group),
        main = 'with RUVs')

normCountData <- normCounts(set_s)
selectedGenes <- names(sort(apply(normCountData, 1, var), 
                            decreasing = TRUE))[1:500]
pheatmap(normCountData[selectedGenes,], 
         annotation_col = colData, 
         show_rownames = FALSE, 
         cutree_cols = 2,
         scale = 'row')

##
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ group)
dds <- dds[rowSums(counts(dds)) > 10, ]
colData(dds) <- cbind(colData(dds),
                      pData(set_s)[rownames(colData(dds)),
                                   grep('W_[0-9]',
                                   colnames(pData(set_s)))])

design(dds) <- ~ W_1 + W_2 + W_3 + group 
dds <- DESeq(dds)
# extract deseq results 
res <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
res <- res[order(res$padj),]

normalized_counts <- counts(dds, normalized = T)
de_genes <- rownames(res)[1:500]
pheatmap(normalized_counts[de_genes,],
         scale = 'row', 
         annotation_col = colData, 
         show_rownames = F, cutree_cols = 2)
