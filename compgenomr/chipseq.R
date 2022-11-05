data_path = system.file('extdata/chip-seq',package='compGenomRData')
chip_files = list.files(data_path, full.names=TRUE)
chip_files

## get chr length
library(GenomeInfoDb)
hg_chrs <- getChromInfoFromUCSC('hg38')
hg_chrs <- subset(hg_chrs, grepl('chr21$', chrom))
seqlengths <- with(hg_chrs, setNames(size, chrom))
seqlengths

library(GenomicRanges)
tilling_window = tileGenome(seqlengths, tilewidth = 1000)
tilling_window = unlist(tilling_window)
tilling_window

library(GenomicAlignments)
bam_files <- list.files(data_path, full.names = T, pattern = 'bam$')
so = summarizeOverlaps(tilling_window, bam_files)
counts = assays(so)[[1]]

cpm = t(t(counts) * (10**6/colSums(counts)))
cpm = cpm[rowSums(cpm) > 0,]
colnames(cpm) = sub('.chr21.bam', '', colnames(cpm))
colnames(cpm) = sub('GM12878_hg38_', '', colnames(cpm))
colnames(cpm)

correlation_mat = cor(cpm, method = 'pearson')
library(ComplexHeatmap)
library(circlize)

heatmap_col = circlize::colorRamp2(breaks = c(-1, 0, 1), 
                                   colors = c('blue', 'white', 'red'))
Heatmap(correlation_mat, heatmap_col)

############
## + and - strand cross-correlations
# load reads
bam_files = list.files(path = data_path, full.names = T, pattern = '.bam$')
bam_files
chip_file = bam_files[7]
reads = readGAlignments(chip_file)
reads = granges(reads)
reads = resize(reads, width=1, fix='start')
reads = keepSeqlevels(reads, 'chr21', pruning.mode = 'coarse')
reads = split(reads, strand(reads))
reads

cov = lapply(reads, function(x){
  coverage(x, width = seqlengths)[[1]] > 0
})

cov <- lapply(cov, as.vector)

wsize = 1:500
jaccard = function(x, y) sum((x & y)) / sum((x|y))

cc = shiftApply(
  SHIFT = wsize,
  X = cov[['+']],
  Y = cov[['-']],
  FUN = jaccard
)

cc = data.frame(fragment_size = wsize, cross_correlation = cc)

library(ggplot2)
ggplot(data = cc, aes(fragment_size, cross_correlation)) +
  geom_point() +
  geom_vline(xintercept = which.max(cc$cross_correlation), 
             size=2, color='red', linetype=2) +
  theme_bw() +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Shift in base pairs') +
  ylab('Jaccard similarity') 

### GC bias identification
library(GenomeInfoDb)
library(GenomicRanges)

hg_chrs <- getChromInfoFromUCSC('hg38')
hg_chrs <- subset(hg_chrs, grepl('chr21$', chrom))
seqlenghts <- with(hg_chrs, setNames(size, chrom))

tilling_window <- unlist(tileGenome(seqlengths, tilewidth = 1000))     
tilling_window

library(BSgenome.Hsapiens.UCSC.hg38)
seq = getSeq(BSgenome.Hsapiens.UCSC.hg38, tilling_window)

nuc = oligonucleotideFrequency(seq, width = 2)
nuc = as.data.frame(nuc)
nuc = round(nuc/1000, 3)

so = summarizeOverlaps(tilling_window, bam_files)
counts = assays(so)[[1]]

cpm = t(t(counts) * (10**6/colSums(counts)))
cpm_log = log10(cpm + 1)

gc = cbind(data.frame(cpm_log), GC = nuc['GC'])
gc

ggplot(gc, aes(GC, GM12878_hg38_CTCF_r1.chr21.bam)) + 
  theme_bw() + geom_point(alpha=0.2)

################ Check the annotations
library(AnnotationHub)
hub = AnnotationHub()

AnnotationHub::query(
  x = hub,
  pattern = c('Ensembl', 'Homo', 'GRCh38', 'chr', 'gtf')
)

gtf = hub[['AH61126']]

ensemble_seqlevels = seqlevels(gtf)
ucsc_seqlevels = paste0('chr', ensemble_seqlevels)
seqlevels(gtf, pruning.mode = 'coarse') = ucsc_seqlevels

gtf = gtf[seqnames(gtf)=='chr21']

annotation_list = GRangesList(
  tss = promoters(x = subset(gtf, type=='gene'),
                  upstream = 1000,
                  downstream = 1000),
  exon= subset(gtf, type=='exon'),
  intron = subset(gtf, type=='gene')
)

bam = readGAlignments(bam_files[1])
results = as.data.frame(
  findOverlaps(bam, annotation_list)
)
head(results)
annotation_name = names(annotation_list)[results$subjectHits]
results$annotation = annotation_name
head(results)
results = results[order(results$subjectHits),]
results = subset(results, !duplicated(queryHits))
library(dplyr)
results = group_by(.data = results, annotation)
results = summarise(.data = results, counts = length(annotation))

results = rbind(results, 
                data.frame(annotation = 'intergenic',
                           counts = length(bam) - sum(results$counts))
)
results$frequency = with(results, round(counts/sum(counts), 2))
results
