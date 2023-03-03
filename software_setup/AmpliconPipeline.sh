conda activate your env
#For Jumana and Reagan
#follow analysis in : https://www.nature.com/articles/s41598-020-71117-4


#1. Download reads to your folder
export BioProject=PRJNA563869
esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv

#Make a test file for just this session:
#You'll want to modify the fastq dump file path for what you set up with Kristina
>gets10
for A in `head $BioProject.SRR`;do
echo "/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 $A" >> gets10; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets10 -n gets10 -a IBN21018 -e cbscott@utexas.edu -t 02:00:00 -w 10 -N 1 -q development
sbatch gets10.slurm

#Make sure it is running
squeue -u cbscott #your username here

#While you wait:
#need to make a new environment for qiime2, seems to have lots of conflicts
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-linux-conda.yml
conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-linux-conda.yml #this is lengthy!!!
conda activate qiime2-2023.2

#We're going to try to use this tutorial: https://docs.qiime2.org/2023.2/
#This data should already be trimmed with the primers removed - let's assume that and keep going
https://docs.qiime2.org/2023.2/tutorials/overview/ #denoising and clustering is where we'll start
https://docs.qiime2.org/2023.2/tutorials/otu-clustering/

#1. Import our data into a qiime artifact rather than fastq files
qiime tools import --input-path SRR10066683.fastq --output-path seqs.qza --type 'SampleData[SequencesWithQuality]'



###after this session do the full batch, 166 files total
>gets_full
for A in `cat $BioProject.SRR`;do
echo "/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 $A" >> gets; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets -n gets -a IBN21018 -e cbscott@utexas.edu -t 08:00:00 -w 24 -N 1 -q normal
