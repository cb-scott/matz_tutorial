#Only need to run these on the first time 
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")

library("dada2")

path <- "~/Downloads/pairedfastq" #be sure to unzip the file 
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
