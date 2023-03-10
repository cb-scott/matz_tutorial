conda activate yourenv
########## Did you sucessfully map reads? ############
#be sure to change things throughout with *your* filepaths! Your files
#will not be in the same location as mine!

#if not, start here:

#### This code creates a test set of fastqs to try mapping
#these will only have 25k reads each and we'll only do 10 files

ls *.fastq | head > sample_fastqs
cat sample_fastqs #should have 10 entries


for file in `cat sample_fastqs`; do
  head -100000 $file > ${file/.fastq}.short.fastq; done

###### Go ahead and map these fastqs quickly to the genome
export GEN=$SCRATCH/your_working_directory/GCA_013753865.1_Amil_v2.1_genomic.fna

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



############START HERE IF YOU HAVE BAM FILES! #########
############# Let's download angsd ################
conda install angsd #this can be touchy - let's see if it works on the first try.

export GEN=$SCRATCH/your_working_directory/GCA_013753865.1_Amil_v2.1_genomic.fna
ls *.sorted.bam > test.sortedbamlist

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -doHWE 1 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
echo "angsd -b test.sortedbamlist -GL 1 $FILTERS $TODO -P 1 -out angsdtest" > angsd_test
python3 $HOME/bin/ls6_launcher_creator.py -j angsd_test -n angsd_test -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q normal

#expected result:
#many files that look lik angsdtest.xxx.gz
#we're looking specifically for a .ibsMat file
