##Load libraries and dependencies
library("Biobase")
library("Biostrings")
library("rtracklayer")
library("GenomeInfoDb")
library("GenomicRanges")
library("Rsamtools")
library("IRanges")
library("GenomicAlignments")
library("ShortRead")
library("chipseq")
library("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome")
library("Gviz")
library("biomaRt")

##---------IMPORT/VIEW BEDFILES------------ 

##set directory to object
dataDirectory = "C:/Users/Owner/Desktop/Bioinformatics_with_R/Datasets/bedfiles"
dataDirectory 

##set BED files as objects (default asRangedData=FALSE)
input = import.bed(file.path(dataDirectory, 'ES_input_filtered_ucsc_chr6.bed'))
rep1 = import.bed(file.path(dataDirectory,'H3K27ac_rep1_filtered_ucsc_chr6.bed'))
rep2 = import.bed(file.path(dataDirectory, 'H3K27ac_rep2_filtered_ucsc_chr6.bed'))


##view data
input
rep1
rep2

##view lengths of each set
length(input)
length(rep1)
length(rep2)

##Questions: what is the length of rep1? What can you infer from the data based on the outputs of input, rep1, and rep2? 

##---------------PREP CHiPseq and controls---------------

##Prepare a single function named prepareChIPseq to perform the read extension
#Estimate mean read length, assign to frag.len
#Extend reads to inferred read length using resize
prepareChIPseq = function(reads){
  frag.len = median(estimate.mean.fraglen(reads))
  cat(paste0('Median fragment size for this library is ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  return(trim(reads.extended))
}


##apply it to the input and CHiP-seq samples
input = prepareChIPseq(input)
rep1 = prepareChIPseq(rep1)
rep1 = prepareChIPseq(rep2)

##Compare with the original input, rep1, and rep2 interval sizes
##Questions: what are the median lengths of the input, rep1, and rep2? 

##------------Bin the CHiP-seq and control---------------
##assign BSgenome annotated mouse mm9 genome to object genome
genome = BSgenome.Mmusculus.UCSC.mm9

##view genome info
genome

##obtain sequence information
seqinfo(genome)


##assign sequence information of BSgenome to object si
si = seqinfo(genome)

##add prefix 'chr' for chromosome names, view chromosome sizes 
si = si[paste0('chr', c(1:19, 'X', 'Y'))]

##view seqinfo object
si

##use tileGenome function to generate a GRanges object in 200bp intervals
binsize = 200
bins = tileGenome(si['chr6'], tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
bins

##count how many reads fall into each bin 
BinCHiPSeq = function(reads, bins) {
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}

##apply it to input, rep1, rep2 
input.200bins = BinCHiPSeq(input, bins)
rep1.200bins = BinCHiPSeq(rep1, bins)
rep2.200bins = BinCHiPSeq(rep2, bins)

rep1.200bins


##plot the coverage for 1000 bins, starting from bin 200,000
plot(200000:201000, rep2.200bins$score[200000:201000], xlab="chr6", ylab="counts per bin")

##export binned data into bedGraph files on local host 
export(input.200bins, con = 'input_chr6.bedGraph', format = "bedGraph")
export(rep1.200bins, con = 'H3K27ac_rep1_chr6.bedGraph', format = "bedGraph")
export(rep2.200bins, con = 'H3K27ac_rep2_chr6.bedGraph', format = "bedGraph")


##-------------Visualization of CHiP-seq data-----------------

##obtain object bm for mm9 from the ensembl 2012 archive 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
               dataset = "mmusculus_gene_ensembl", 
               host = "http://may2012.archive.ensembl.org")
fm = Gviz:::.getBMFeatureMap()
fm["symbol"] = "external_gene_id"

##view specifically chromosome 6 starting at designated start/end position, this region encodes a highly ES cell specific Nanog gene
##isolate gene models for this interval, save to data directory as object bm
bm = BiomartGeneRegionTrack(chromosome = 'chr6', genome = "mm9", 
                            start=122530000, end=122900000, 
                            biomart = mart, filters = list("with_ox_refseq_mrna"=TRUE),
                            size=4, name = "RefSeq", utr5="red3", utr3="red3",
                            protein_coding="black", col.line=NULL, cex=7,
                            collapseTranscripts="longest",
                            featureMap = fm)
##view bm 
bm

##assign genomeaxistrack to object AT 
AT = GenomeAxisTrack()

##plot results using plotTracks
plotTracks(c(bm, AT), 
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window = "auto",
           cex.title=1, fontsize=10)

##What we have here (you can view in plots tab) is the reference mouse genome mm9 at chromosome 6 for our given genes. 
##Next we need to add our sample tracks for comparison. 

input.track = DataTrack(input.200bins, strand = "*", genome = "mm9", col.histogram='gray',
                        fill.histogram='black', name = "input", col.axis="black",
                        cex.axis=0.4, ylim=c(0,150))
rep1.track = DataTrack(rep1.200bins, strand = "*", genome = "mm9", col.histogram='steelblue',
                        fill.histogram='black', name = "Rep. 1", col.axis="steelblue",
                        cex.axis=0.4, ylim=c(0,150))

rep2.track = DataTrack(rep2.200bins, strand = "*", genome = "mm9", col.histogram='steelblue',
                       fill.histogram='black', name = "Rep. 2", col.axis="steelblue",
                       cex.axis=0.4, ylim=c(0,150))
##plot sample tracks
plotTracks(c(input.track, rep1.track, rep2.track, bm, AT), 
           from = 122530000, to = 122900000,
           transcriptAnnotation="symbol", window="auto",
           type="histogram", cex.title=0.7, fontsize = 10)

##View the plots that have been generated. What inferences can you make about the gene expression of our sample datasets? 


##-----------------PEAK CALLING----------------

##import files that have already been run through MACS, with isolated peaks 
peaks.rep1 = import.bed(file.path(dataDirectory,'Rep1_peaks_ucsc_chr6.bed'))
peaks.rep2 = import.bed(file.path(dataDirectory, 'Rep2_peaks_ucsc_chr6.bed'))

##display peaks in our browsers as blue boxes 
peaks1.track = AnnotationTrack(peaks.rep1, genome = "mm9", name = 'Peaks Rep. 1',
                               chromosome = 'chr6', shape='box', fill='blue3', size=2)

peaks2.track = AnnotationTrack(peaks.rep2, genome = "mm9", name = 'Peaks Rep. 2',
                               chromosome = 'chr6', shape='box', fill='blue3', size=2)

##visualize the Nanog locus 
plotTracks(c(input.track, rep1.track, peaks1.track, rep2.track, peaks2.track, bm, AT),
           from = 122530000, to = 122900000,
           transcriptAnnotation = "symbol", window = "auto", 
           type = "histogram", cex.title=0.7, fontsize=10)      



##This is the end of this workshop, although we didn't cover much today I hope you enjoyed it! 
##Pretty much everything here was pulled from 
##https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-5-chipseq/Epigenetics.html
##If you would like to know how to do more downstream analyses, such as peak calling and further data visualizations going through the entire page is a great start


