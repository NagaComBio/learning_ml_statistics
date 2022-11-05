library(compGenomRData)
library(tidyverse)

counts_file <-  system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                            package='compGenomRData')
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package='compGenomRData')

counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
head(counts)
dim(counts)

## calculate TPM
genelength <- counts[,'width']
counts <- counts[,1:10]
counts_rwgl <- apply(counts, 2, function(x) x / (genelength/1000))
counts_tmp <- apply(counts_rwgl, 2, function(x) (x / sum(x)) * 10**6)
counts_tmp[,1:10]

colSums(counts_tmp)

### high varianbles
library(pheatmap)
var_genes <- apply(counts_tmp, 1, var) 
high_var_genes <- names(sort(var_genes, decreasing = T)[1:100])

pheatmap(counts_tmp[high_var_genes,], scale = 'row')

## Sample correlations
cor(counts_tmp)
corrplot::corrplot(cor(counts_tmp), method = 'ellipse', 
                   type = 'upper', 
                   hclust.method = 'average')
## highest exp
high_expre_genes <- names(sort(rowSums(counts_tmp), decreasing = T))[1:100]
pheatmap(counts_tmp[high_expre_genes,], scale = 'row')
library(ggfortify)
df <- as.data.frame(t(counts_tmp[high_expre_genes,]))
pca <- prcomp(df, scale. = T)
df %>% rownames_to_column('sample_name') %>%
  mutate(sample_name = factor(sample_name)) %>%
  mutate(group= ifelse(grepl('CASE', sample_name), 'Case', 'Control')) -> df

autoplot(pca, data = df, colour='group')

t(counts_tmp[high_expre_genes,])

library(Rtsne)
tsne <- Rtsne(t(counts_tmp[high_expre_genes,]), perplexity = 2)
tnse<-as.data.frame(tsne$Y)
tnse$sample_nsme <- df$sample_name
tnse$group <- df$group
ggplot2::ggplot(tnse, aes(V1, V2, color=group)) + geom_point()