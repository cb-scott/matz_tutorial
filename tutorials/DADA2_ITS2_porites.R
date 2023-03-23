library(dada2)
library(ShortRead)
library(Biostrings)

#chekc for primer seqs 
path <- "data/POR_ITS2/"
list.files(path)

fnFs <- sort(list.files(path, pattern = "1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "2.fastq.gz", full.names = TRUE))


#used cudadapt to remove primers already, so skip over

# Forward and reverse fastq filenames have the format:
cutFs <-fnFs
cutRs <- fnRs

#ITS2
FWD="GAATTGCAGAACTCCGTGAACC"
REV="CGGGTTCWCTTGTYTGACTTCATGC"

#Illumina
#FWD="CTACACGACGCTCTTCCGATCT"
#REV="CAGACGTGTGCTCTTCCGATCT"


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][2]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])


fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

DA#check for primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#primers have been successfully trimmed....and illumina adapters....



path.cut = "data/POR_ITS2/filtN/"
cutFs <- sort(list.files(path.cut, pattern = "1.fastq.gz", full.names = T))
cutRs <- sort(list.files(path.cut, pattern = "2.fastq.gz", full.names = T))


filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
#takes forever!
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 5), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set mult

path.filts = "data/POR_ITS2/filtN/filtered/"
filtFs <- sort(list.files(path.filts, pattern = "1.fastq.gz", full.names = T))
filtRs <- sort(list.files(path.filts, pattern = "2.fastq.gz", full.names = T))


#these all take forever - definitely needs to be on TACC?
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = T)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
sample.names <- unname(sapply(filtFs, get.sample.name)) 

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Considerations for your own data: Most of your reads should successfully merge.
#If that is not the case upstream parameters may need to be revisited: 
#Did you trim away the overlap between your reads? Potentially!!!
#Might have primers at start and end of read inflating sequences... may want to pass cutadapt again?
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate = T) #just concatenate doesn't consider paired end 

#I don't think my reads overlap - after trimming primers read length only 126, target region is longer
#for more detail:
#https://github.com/benjjneb/dada2/issues/790
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 2)

seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#grab a database from Danielle Claar
#https://github.com/baumlab/Claar_etal_2020_SciRep/blob/master/ITS2db_trimmed_derep_dada.fasta

its2.ref <- "data/ITS2db_trimmed_derep_dada.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, its2.ref, multithread = TRUE, tryRC = TRUE)


taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tp <- data.frame(taxa.print)
tp$Genus
unique(tp$Species)

###########################################
## what if we don't merge? just classify F and R
# 
# seqtab <- makeSequenceTable(dadaRs)
# table(nchar(getSequences(seqtab)))
# 
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)
# 
# #grab a database from Danielle Claar
# #https://github.com/baumlab/Claar_etal_2020_SciRep/blob/master/ITS2db_trimmed_derep_dada.fasta
# 
# its2.ref <- "data/ITS2db_trimmed_derep_dada.fasta"  # CHANGE ME to location on your machine
# taxa <- assignTaxonomy(seqtab.nochim, its2.ref, multithread = TRUE, tryRC = TRUE)
# 
# 
# taxa.print <- taxa  # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)
# 
# tp <- data.frame(taxa.print)
# tp$Genus
# table(tp$Species)
# table(tp$Genus)
# 
# #we get better resolution?!
# 


##########################################
#### GO WITH THE MERGED PAIRS FOR NOW#####
##########################################
#see: https://www.biostars.org/p/9498358/

#BiocManager::install("phyloseq")
library(phyloseq)
