#GET DESEQ COUNT -> PCA TUTORIAL
#The first time you run this, download DESeq2
#DESeq2 is a software that was originally built for gene expression data
#that normalizes the read counts per gene and normalizes based
#on sequencing effort. We have a super similar problem here - we have different 
#read counts per ASV, so we want to normalize them, calculate distances
#then create a PCA.

#Phyloseq is still great for looking at abundances of different taxa, but
#we'll just make the PCA by hand here for simplicity!

#install DESeq2 the first time you run this
BiocManager::install("DESeq2")

#after it's installed, start here
library(DESeq2)
library(phyloseq)
library(vegan)
library(pheatmap)
library(tidyverse)

#SAVE your phyloseq object as an RData object from your last script. 
#I think we already did this.

#read it in here. Substitute my path with yours. This object should have otu_tables, metadata, etc
ps = readRDS("data/ITS2_ps.rds")

#convert it to a deseq2 object
its2.dds = phyloseq_to_deseq2(ps, ~ 1)

#Do some quality control
startdim = dim(counts(its2.dds))
keeprows <- rowSums(counts(its2.dds)) >= 10 #remove any ASV's that don't have at least 10 reads mapping to them
keepcols <- colSums(counts(its2.dds)) >=50 #trash samples with fewer than 50 total read count
its2.dds <- its2.dds[keeprows,keepcols]
enddim = dim(counts(its2.dds))

#how many samples did we loose?
#Stop here if this is more than 20-25%!!! Talk to me about better filtering
print(paste("Starting ASVs=", startdim[1], "Ending ASVS=", enddim[1]))
print(paste("Starting Samples=", startdim[2], "Ending Samples=", enddim[2]))

#do a deSeq Transformation to estimate sample variance
its2.dds <- estimateSizeFactors(its2.dds, type = 'poscounts')

#Get Variance Stabilized Counts for your ASVs. These counts control for
#sampling effort and the variance between your samples. You may get
#negative numbers and decimals here - that is ok!
vsd <- varianceStabilizingTransformation(its2.dds, blind = T)
count.tab = t(assay(vsd))
#count.tab[count.tab <0] = 0
#Calculate the distance between your samples. 
sampleDists <- vegdist(count.tab, method = "manhattan")

#create a heatmap of your samples - this distance is not constrained between 0 and 1
pheatmap(sampleDists)

#Now make a pca!
#Capscale is the PCA funciton from Vegan
cap = capscale(sampleDists~1)
plot(cap)

#take a look at it with respect to your metadata 
#you may have already read this in up top, if so, skip likes 65 and 66!
meta = sample_data(ps)
meta$ID = rownames(meta) 

pc_scores = data.frame(scores(cap)$sites)
pc_scores$ID = rownames(pc_scores)

#combine your pca scores with your metadata 
pc_meta = left_join(pc_scores, meta, by=c("ID" = "ID"))
head(pc_meta) #did it work?

#now plot a nice pca using ggplot
#Plot this in several different ways looking at your variables! Simply change
#the col = Age argument to whatever you want to plot from your data. 

#what are your metadata columns?
colnames(meta)
#instead of age and site below, replace with column names from your data 

pc_meta %>% ggplot(aes(x = MDS1, y = MDS2, col = Age)) + geom_point(size = 3, alpha = .8) + 
  theme_classic() + coord_fixed() + geom_vline(xintercept = 0, lty = 'dashed', col = 'grey') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey')


pc_meta %>% ggplot(aes(x = MDS1, y = MDS2, col = Site)) + geom_point(size = 3, alpha = .8) + 
  theme_classic() + coord_fixed() + geom_vline(xintercept = 0, lty = 'dashed', col = 'grey') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey')

#change shape of points
pc_meta %>% ggplot(aes(x = MDS1, y = MDS2, col = Site, shape = Age)) + geom_point(size = 3, alpha = .8) + 
  theme_classic() + coord_fixed() + geom_vline(xintercept = 0, lty = 'dashed', col = 'grey') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey')


