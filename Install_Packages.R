##install packages not associated with bioconductor 
install.packages("stringi")

##install bioconductor base and dependencies
source("https://bioconductor.org/biocLite.R")
##update bioconductor base
biocLite()

##Additional packages in Bioconductor: https://www.bioconductor.org/packages/release/BiocViews.html#___Software
BiocInstaller::biocLite("Rsamtools")
BiocInstaller::biocLite("GenomicRanges")
BiocInstaller::biocLite("GenomeInfoDb")
#BiocInstaller::biocLite("GEOquery")
#BiocInstaller::biocLite("seqinr")
#BiocInstaller::biocLite("SRAdb")
BiocInstaller::biocLite("ShortRead")
BiocInstaller::biocLite("rtracklayer")
#BiocInstaller::biocLite("VariantAnnotation")
#BiocInstaller::biocLite("GenomicFeatures")
BiocInstaller::biocLite("Biostrings")
BiocInstaller::biocLite("IRanges")
BiocInstaller::biocLite("BSgenome")
BiocInstaller::biocLite("chipseq")
BiocInstaller::biocLite("SummarizedExperiment")
BiocInstaller::biocLite("BSgenome.Mmusculus.UCSC.mm9")
BiocInstaller::biocLite("Gviz")
BiocInstaller::biocLite("biomaRt")
BiocInstaller::biocLite("seqinfo")
