conda activate your_env
####Grab and build a reference genome
#This is a bit separate from Kristina's project, so let's put it all in a new directory
#Let's do Acropora Millepora
conda install ncbi-datasets-cli
#(for Carly - conda activate metagenomics)

datasets download genome accession GCA_013753865.1 # can run this on login line, should take a couple minutes

#unzip the folder
unzip ncbi_dataset.zip
cd ncbi_dataset
cd data
cd GCA_013753865.1/

#look at what your genome looks like
head GCA_013753865.1_Amil_v2.1_genomic.fna
mv GCA_013753865.1_Amil_v2.1_genomic.fna $SCRATCH/your_working_directory
#Now let's build this! We'll need bowtie2 - which I think you have installed?
#Let me know if this part gives you an error
export GEN=$SCRATCH/qiimetest/GCA_013753865.1_Amil_v2.1_genomic.fna #THE FULL PATH WHERE YOUR FILE IS LOCATED!
echo "bowtie2-build $GEN $GEN" >btb
ls6_launcher_creator.py -j btb -n btb -a IBN21018 -e YOUREMAIL -t 02:00:00 -N 1 -w 1 -q normal
sbatch btb.slurm
#check that it batched!
squeue -u cbscott # your username here

########Grab reads while this is downloading##################
#use bioproject PRJNA434194

#give this a try on your own!

export BioProject=PRJNA434194
esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv

#Make a test file for just this session:
#You'll want to modify the fastq dump file path for what you set up with Kristina
>gets10
for A in `head $BioProject.SRR`;do
echo "/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 $A" >> gets10; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets10 -n gets10 -a IBN21018 -e YOUREMAIL -t 00:20:00 -w 10 -N 1 -q development
sbatch gets10.slurm
squeue -u cbscott #check if it's running


###after this session do the full batch - way more files ######
>gets_full
for A in `cat $BioProject.SRR`;do
echo "/work/06909/cbscott/ls6/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump-orig.3.0.0 $A" >> gets; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets -n gets -a IBN21018 -e YOUREMAIL -t 08:00:00 -w 24 -N 1 -q normal

######Map reads, exciting!######
>map_amil
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN -U $F -S ${F/.trim}.master.local.sam && \
samtools sort -O bam -o ${F/.trim}.master.local.sorted.bam ${F/.trim}.master.local.sam && samtools index -c ${F/.trim}.master.local.sorted.bam" >> map_amil; done
python3 $HOME/bin/ls6_launcher_creator.py -j map_amil -n map_dolphin -a IBN21018 -e YOUREMAIL -t 01:00:00 -N 1 -w 10 -q development
sbatch map_amil.slurm
