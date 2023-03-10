################# INSTALLATION INSTRUCTIONS #######################
#need to make a new environment for qiime2, seems to have lots of conflicts
#Only need to do this the first time
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-linux-conda.yml
conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-linux-conda.yml #this is lengthy!!!
conda activate qiime2-2023.2 #check it
######################################################################
#activate your old environment to do this processing
conda activate your_old_environment_from_our_setup
cds #change into SCRATCH
mkdir microbiome_test #or whatever you want to call it - make a clean folder for this.
cd microbiome_test
mkdir pairedfastq #call it this so the rest of the code works
cd pairedfastq
#1. Download reads to your folder, run in 'pairedfastq'
export BioProject=PRJNA563869
esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv

#Make a test file for just this session:
#You'll want to modify the fastq dump file path for what you set up with Kristina
#/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 carly's path
>gets10
for A in `head $BioProject.SRR`;do
echo "fastq-dump-orig.3.0.1 --split-files $A" >> gets10; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets10 -n gets10 -a IBN21018 -e cbscott@utexas.edu -t 00:10:00 -w 10 -N 1 -q development
sbatch gets10.slurm

#https://rstudio-pubs-static.s3.amazonaws.com/794628_9729e7e28b0049ab98b61adedbdda2a5.html
conda activate qiime2-2023.2

rename -v '_' '_00_L001_' * /
rename -v '.fastq' '_001.fastq' * /
rename -v '_1' '_R1' * /
rename -v '_2' '_R2' *
gzip *.fastq
#this works!
cd .. #change directories one up, now should be in $SCRATCH/microbiome_test (or whatever you called it)
#it's a little slow so,
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path pairedfastq --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

### Feel free to try to work ahead based on the "otu-clustering" tutorial below if you feel like it.
echo "qiime vsearch dereplicate-sequences --i-sequences demux-paired-end.qza --o-dereplicated-table table.qza --o-dereplicated-sequences rep-seqs.qza" > derep
python3 $HOME/bin/ls6_launcher_creator.py -j derep -n derep -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -w 1 -N 1 -q development


#We're going to try to use this tutorial: https://docs.qiime2.org/2023.2/
#This data should already be trimmed with the primers removed - let's assume that and keep going
https://docs.qiime2.org/2023.2/tutorials/overview/ #denoising and clustering is where we'll start
https://docs.qiime2.org/2023.2/tutorials/otu-clustering/


###DONT do this yet - let's get it workign first.
>gets_full
for A in `cat $BioProject.SRR`;do
echo "/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 $A" >> gets; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets -n gets -a IBN21018 -e cbscott@utexas.edu -t 08:00:00 -w 24 -N 1 -q normal
