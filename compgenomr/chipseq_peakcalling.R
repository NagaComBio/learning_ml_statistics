data_path = system.file('extdata/chip-seq',package='compGenomRData')
# set names for chip-seq bigWig files
chip_files = list(
  H3K4me3  = 'GM12878_hg38_H3K4me3.chr21.bw',
  
  H3K36me3 = 'GM12878_hg38_H3K36me3.chr21.bw',
  
  POL2     = 'GM12878_hg38_POLR2A.chr21.bw'
)
# get full paths to the files
chip_files = lapply(chip_files, function(x){
  file.path(data_path, x)
})

chip_files
## Load bw profiles
library(rtracklayer)
chip_profiles <- lapply(chip_files, rtracklayer::import.bw)
chip_profiles

## load annotation
library(AnnotationHub)
hub = AnnotationHub()
gtf = hub[['AH61126']]

seqlevels(gtf, pruning.mode='coarse') = '21'
ensembl_seqlevels = seqlevels(gtf)
ucsc_seqlevels = paste0('chr', ensembl_seqlevels)
seqlevels(gtf, pruning.mode='coarse') = ucsc_seqlevels 

## genomic ranges to transcript database
library(GenomicFeatures)
library(Gviz)
txdb <- makeTxDbFromGRanges(gtf)
gene_track = GeneRegionTrack(txdb, chr='chr21', genome = 'hg38')

hg_chrs <- getChromInfoFromUCSC('hg38')
hg_chrs = subset(hg_chrs, (grepl('chr21$', chrom)))
seqlengths = with(hg_chrs, setNames(size, chrom))

chr_track = IdeogramTrack(
  chromosome = 'chr21',
  genome = 'hg38'
)

axis = GeneRegionTrack(
  range = GRanges('chr21', IRanges(1, width = seqlengths))
)

data_tracks = lapply(names(chip_profiles), function(exp_name){
  
  # chip_profiles[[exp_name]] - selects the 
  # proper experiment using the exp_name
  DataTrack(
    range = chip_profiles[[exp_name]],   
    name  = exp_name,  
    
    # type of the track
    type  = 'h', 
    
    # line width parameter
    lwd   = 5
  )
})


# select the start coordinate for the URB1 gene
start = min(start(subset(gtf, gene_name == 'URB1')))

# select the end coordinate for the URB1 gene
end   = max(end(subset(gtf, gene_name == 'URB1')))

# plot the signal profiles around the URB1 gene
plotTracks(
  trackList = c(chr_track, axis, gene_track, data_tracks),
  
  # relative track sizes
  sizes     = c(1,1,1,1,1,1), 
  
  # background color
  background.title     = "black",
  
  # controls visualization of gene sets
  collapseTranscripts  = "longest", 
  transcriptAnnotation = "symbol",
  
  # coordinates to visualize 
  from = start - 5000,
  to   = end   + 5000
)

gene_track

###########################
## Narrow peak calling

library(normR)
# full path to the ChIP data file
chip_file    = file.path(data_path, 'GM12878_hg38_CTCF_r1.chr21.bam')

# full path to the Control data file
control_file = file.path(data_path, 'GM12878_hg38_Input_r5.chr21.bam')

# as previously done, we calculate the cpm for each experiment
library(GenomicRanges)
library(GenomicAlignments)

# select the chromosome
hg_chrs = getChromInfoFromUCSC('hg38') 
hg_chrs = subset(hg_chrs, grepl('chr21$',chrom))

seqlengths = with(hg_chrs, setNames(size, chrom))

# define the windows
tilling_window = unlist(tileGenome(seqlengths, tilewidth=1000))

# count the reads
counts         = summarizeOverlaps(
  features = tilling_window, 
  reads    = c(chip_file, control_file)
)

# normalize read counts
counts         = assays(counts)[[1]]
cpm = t(t(counts)*(1000000/colSums(counts)))

library(ggplot2)
# convert the matrix into a data.frame for ggplot
cpm = data.frame(cpm)
ggplot(
  data = cpm, 
  aes(
    x = GM12878_hg38_Input_r5.chr21.bam, 
    y = GM12878_hg38_CTCF_r1.chr21.bam)
) +
  geom_point() +
  geom_abline(slope = 1) +
  theme_bw() +
  theme_bw() +
  scale_fill_brewer(palette='Set2') +
  theme(
    axis.text   = element_text(size=10, face='bold'),
    axis.title  = element_text(size=14,face="bold"),
    plot.title  = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Input CPM') +
  ylab('CTCF CPM') +
  ggtitle('ChIP versus Input')

library(normr)
# peak calling using chip and control
ctcf_fit = enrichR(
  
  # ChIP file
  treatment = chip_file,
  
  # control file
  control   = control_file,
  
  # genome version
  genome    = "hg38",
  
  # print intermediary steps during the analysis
  verbose   = FALSE)

summary(ctcf_fit)
# extracts the ranges
ctcf_peaks = getRanges(ctcf_fit)

# annotates the ranges with the supporting p value
ctcf_peaks$qvalue     = getQvalues(ctcf_fit)

# annotates the ranges with the calculated enrichment
ctcf_peaks$enrichment = getEnrichment(ctcf_fit)

# selects the ranges which correspond to the enriched class
ctcf_peaks = subset(ctcf_peaks, !is.na(component))

# filter by a stringent q value threshold
ctcf_peaks = subset(ctcf_peaks, qvalue < 0.01)

# order the peaks based on the q value
ctcf_peaks = ctcf_peaks[order(ctcf_peaks$qvalue)]
ctcf_peaks
seqlevels(ctcf_peaks, pruning.mode='coarse') = 'chr21'

# write the peaks loacations into a txt table
write.table(ctcf_peaks, file.path(data_path, 'CTCF_peaks.txt'), 
            row.names=F, col.names=T, quote=F, sep='\t')

# find enriched tilling windows
enriched_regions = countOverlaps(tilling_window, ctcf_peaks) > 0
library(ggplot2)
cpm$enriched_regions = enriched_regions 

ggplot(
  data = cpm, 
  aes(
    x = GM12878_hg38_Input_r5.chr21.bam, 
    y = GM12878_hg38_CTCF_r1.chr21.bam, 
    color = enriched_regions
  )) +
  geom_point() +
  geom_abline(slope = 1) +
  theme_bw() +
  scale_fill_brewer(palette='Set2') +
  theme(
    axis.text   = element_text(size=10, face='bold'),
    axis.title  = element_text(size=14,face="bold"),
    plot.title  = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Input CPM') +
  ylab('CTCF CPM') +
  ggtitle('ChIP versus Input') +
  scale_color_manual(values=c('gray','red'))
############################
## Broad peak calling

# fetch the ChIP-file for H3K36me3
chip_file    = file.path(data_path, 'GM12878_hg38_H3K36me3.chr21.bam')

# fetch the corresponding input file
control_file = file.path(data_path, 'GM12878_hg38_Input_r5.chr21.bam')

library(normr)
countConfiguration = countConfigSingleEnd(binsize = 5000)

# find broad peaks using enrichR
h3k36_fit = enrichR(
  
  # ChIP file
  treatment   = chip_file,
  
  # control file
  control     = control_file,
  
  # genome version
  genome      = "hg38",
  verbose     = FALSE,
  
  # window size for counting
  countConfig = countConfiguration)

summary(h3k36_fit)
h3k36_peaks = getRanges(h3k36_fit)
h3k36_peaks$qvalue = getQvalues(h3k36_fit)
h3k36_peaks$enrichment = getEnrichment(h3k36_fit)
h3k36_peaks<-subset(h3k36_peaks, !is.na(h3k36_peaks$component))
h3k36_peaks <- subset(h3k36_peaks, h3k36_peaks$qvalue < 0.01)
h3k36_peaks <- h3k36_peaks[order(h3k36_peaks$qvalue)]

h3k36_peaks <- reduce(h3k36_peaks)

counts  = summarizeOverlaps(
  features = tilling_window, 
  reads    = c(chip_file, control_file)
)

# normalize read counts
counts         = assays(counts)[[1]]
cpm = t(t(counts)*(1000000/colSums(counts)))
cpm <- as.data.frame(cpm)
ggplot(cpm, 
       aes(GM12878_hg38_Input_r5.chr21.bam,
           GM12878_hg38_H3K36me3.chr21.bam)) + 
  geom_point(alpha=0.2) + 
  geom_abline(slope = 1)

enriched_regions = countOverlaps(tilling_window, h3k36_peaks) > 0
cpm$enriched_regions = enriched_regions 

ggplot(cpm, 
       aes(GM12878_hg38_Input_r5.chr21.bam,
           GM12878_hg38_H3K36me3.chr21.bam,
           color=enriched_regions)) + 
  geom_point(alpha=0.2) + 
  geom_abline(slope = 1)
####################### 
## chip sea peak quality control
h3k36_counts <- as.data.frame(getCounts(h3k36_fit))
head(h3k36_counts)
colnames(h3k36_counts) <- c('Input', 'H3K36me3')
h3k36_counts$qvalue = getQvalues(h3k36_fit)
h3k36_counts$enriched[is.na(h3k36_counts$qvalue)] = "Not Peak"
h3k36_counts$enriched[h3k36_counts$qvalue > 0.05] = "Not Peak"
h3k36_counts$enriched[h3k36_counts$qvalue <= 0.05] = "Peak"
h3k36_counts$qvalue <- NULL

head(h3k36_counts)
h3k36_counts_df = tidyr::pivot_longer(
  data = h3k36_counts,
  cols = -enriched,
  names_to = 'experiment',
  values_to = 'counts'
)
h3k36_counts_df %>% 
  dplyr::group_by(experiment, enriched) %>%
  dplyr::summarize(num_of_reads = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(experiment) %>%
  dplyr::mutate(total = sum(num_of_reads)) %>% 
  dplyr::mutate(percentage = round(num_of_reads/total, 2)) %>% 
  ggplot(aes(experiment, percentage, fill=enriched)) + 
    geom_bar(stat='identity', position = 'dodge')


### ### 
## Check DNA motifs in the peaks
library(MotifDb)

motifs = query(query(MotifDb, 'Hsapiens'), 'CTCF')
ctcf_motif = motifs[[1]]
ctcf_peaks_resized = resize(ctcf_peaks, width=400, fix = 'center')
# load the human genome sequence
library(BSgenome.Hsapiens.UCSC.hg38)
# extract the sequences around the peaks
seq = getSeq(BSgenome.Hsapiens.UCSC.hg38, ctcf_peaks_resized)

library(TFBSTools)
# convert the matrix into a PWM object
ctcf_pwm = PWMatrix(
  ID = 'CTCF', 
  profileMatrix = ctcf_motif
)

sum(apply(ctcf_motif, 2, max))
hits = searchSeq(ctcf_pwm, seq, min.score="80%", strand="*")
# convert the hits into a data.frame
hits = as.data.frame(hits)
motif_hits_df = data.frame(
  peak_order = 1:length(ctcf_peaks)
)
motif_hits_df$contains_motif = motif_hits_df$peak_order %in% hits$seqnames
motif_hits_df = motif_hits_df[order(-motif_hits_df$peak_order),]
head(motif_hits_df)

motif_hits_df$perc_peaks = with(motif_hits_df, 
                                cumsum(contains_motif) / max(peak_order))
motif_hits_df$perc_peaks = round(motif_hits_df$perc_peaks, 2)
# plot the cumulative distribution of motif hit percentages
library(ggplot2)
ggplot(
  motif_hits_df, 
  aes(
    x = peak_order, 
    y = perc_peaks
  )) +
  geom_line(size=2) +
  theme_bw() +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Peak rank') +
  ylab('Percetage of peaks with motif') +
  ggtitle('Percentage of CTCF peaks with the CTCF motif')

### Motif localization
ctcf_peaks_resized = resize(ctcf_peaks, width = 2000, fix='center')
# fetch the sequence
seq = getSeq(BSgenome.Hsapiens.UCSC.hg38,ctcf_peaks_resized)

# convert the motif matrix to PWM, and scan the peaks
ctcf_pwm    = PWMatrix(ID = 'CTCF', profileMatrix = ctcf_motif)
hits = searchSeq(ctcf_pwm, seq, min.score="80%", strand="*")
hits = as.data.frame(hits)
head(hits)
# set the position relative to the start
hits$position = hits$start - 1000 
head(hits)
# plot the motif hits around peaks
ggplot(data=hits, aes(position)) +
  geom_density(size=2) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype=2, color='red', size=2) +
  xlab('Position around the CTCF peaks') +
  ylab('Per position percentage\nof motif occurence') +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5))


#####################
## Annotation of the peaks with gene elements
library(AnnotationHub)
hub = AnnotationHub()
gtf = hub[['AH61126']]

seqlevels(gtf, pruning.mode='coarse') = '21'
seqlevels(gtf, pruning.mode='coarse') = paste0('chr21')

annotation_list = GRangesList(
  tss = promoters(subset(gtf, type=='gene'), 1000, 1000),
  exon = subset(gtf, type=='exon'),
  intron = subset(gtf, type=='gene') 
)
annotation_list

#ctcf_peaks_regions = as.data.frame(findOverlaps(reduce(ctcf_peaks), annotation_list))
ctcf_peaks_regions = as.data.frame(findOverlaps(reduce(h3k36_peaks), annotation_list))

head(ctcf_peaks_regions)

ctcf_peaks_regions$annotation = names(annotation_list)[ctcf_peaks_regions$subjectHits]
lapply(ctcf_peaks_regions, length)
ctcf_peaks_regions = subset(ctcf_peaks_regions, !duplicated(queryHits))

ctcf_peaks_regions %>% 
  dplyr::group_by(annotation) %>% 
  dplyr::summarise(count = length(annotation)) -> ctcf_peaks_regions

ctcf_peaks_regions %>% 
  dplyr::add_row(annotation = 'integenic', 
                 count = length(ctcf_peaks) - sum(ctcf_peaks_regions$count)) -> ctcf_peaks_regions
library(ggplot2)
ctcf_peaks_regions %>%
  dplyr::mutate('peaks' = 'H3H36me3') %>% 
  dplyr::mutate(percentage = count/sum(count)) %>%
  ggplot2::ggplot(aes(peaks, percentage, fill=annotation)) + 
    geom_bar(stat = 'identity') -> p2
