conda activate your_env
cds # work in scratch
mkdir your_working_directory #<- put your files in a new directory

####Grab and build a reference genome
#This is a bit separate from Kristina's project, so let's put it all in a new directory
#Let's do Acropora Millepora
conda install ncbi-datasets-cli # only need to do this the first time.

########DOWLOAD THE GENOME FOR ACROPORA MILLEPORA
datasets download genome accession GCA_013753865.1 # can run this on login line, should take a couple minutes
#unzip the folder
unzip ncbi_dataset.zip
cd ncbi_dataset
cd data
cd GCA_013753865.1/
#look at what your genome looks like
head GCA_013753865.1_Amil_v2.1_genomic.fna
#move the genome to your working directory (you made this up top)
mv GCA_013753865.1_Amil_v2.1_genomic.fna $SCRATCH/your_working_directory
#change into your working direcotry
cd $SCRATCH/your_working_directory
ls #list and make sure that the genome file is in your "present working directory"



#########BUILD YOUR GENOME USING BOWTIE ######################
#Now let's build this! We'll need bowtie2 - which I think you have installed?
#Let me know if this part gives you an error
#THE FULL PATH WHERE YOUR FILE IS LOCATED! You'll have to change this!
export GEN=$SCRATCH/your_working_directory/GCA_013753865.1_Amil_v2.1_genomic.fna
echo "bowtie2-build $GEN $GEN" >btb
ls6_launcher_creator.py -j btb -n btb -a IBN21018 -e YOUREMAIL -t 02:00:00 -N 1 -w 1 -q normal
sbatch btb.slurm
#check that it batched!
squeue -u cbscott # your username here
#wait for job to run. This should take about 10 minutes.
#DO NOT batch an additonal job.
ll -tr #list the contents of your directory, sorted by time
#Expected Output: an additional 4-6 files called something like
GCA_013753865.1_Amil_v2.1_genomic.fna.(somenumber).bt2

########Grab reads while this is downloading##################
#use bioproject PRJNA434194
#RUN this in the same directory as where you put your GENOME
pwd
#this should print $SCRATCH/your_working_directory
export BioProject=PRJNA434194
esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv

#Make a test file for just this session:
>gets10
for A in `head $BioProject.SRR`;do
echo "fastq-dump-orig.3.0.1  $A" >> gets10; done
ls6_launcher_creator.py -j gets10 -n gets10 -a IBN21018 -e YOUREMAIL -t 00:20:00 -w 10 -N 1 -q development
sbatch gets10.slurm
squeue -u cbscott #check if it's running

#wait for your job to finish - <5 minutes
#check to make sure it worked:
#EXPECTED OUTPUT
ll -tr
#you see 10 *.fastq files, created within the last 10 minutes (should have a time signature)

######Map reads, exciting!######
#this aligns your reads to the reference genome
echo $GEN
#this should print the path to your ref genome, if it doesn't go back and check/
#rerun line 31
#Again run this in the same directory that has your *.fastq files and genome.
#####!!!!! You should only have fastq files for this species in this folder!!!!!####
>map_amil
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN -U $F -S ${F/.trim}.master.local.sam && \
samtools sort -O bam -o ${F/.trim}.master.local.sorted.bam ${F/.trim}.master.local.sam && samtools index -c ${F/.trim}.master.local.sorted.bam" >> map_amil; done
ls6_launcher_creator.py -j map_amil -n map_amil -a IBN21018 -e YOUREMAIL -t 01:00:00 -N 1 -w 10 -q development
sbatch map_amil.slurm
#wait for job to finish - this might take a while.
#after job finishes,
ll -tr
#EXPECTED RESULT
#at least 10 files each that end in .bam, .sam, and .csi


###after you get this all working, download all of the files and try mapping ######
>gets_full
for A in `cat $BioProject.SRR`;do
echo "fastq-dump-orig.3.0.1 $A" >> gets; done
python3 $HOME/bin/ls6_launcher_creator.py -j gets -n gets -a IBN21018 -e YOUREMAIL -t 04:00:00 -w 24 -N 1 -q normal
