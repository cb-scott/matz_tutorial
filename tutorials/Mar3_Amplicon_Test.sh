conda activate qiime2-2023.2
#make sure the environment number matches
conda install -c https://packages.qiime2.org/qiime2/2023.2/tested/ -c conda-forge -c bioconda -c defaults q2-fondue



#will need to try to make a sample manifest
echo "sample-id,absolute-filepath" > manifest.csv
>sampleid
for file in *.fastq; do echo ${file/.fastq} >> sampleid; done
readlink -f *.fastq >fullfilepath

paste -d, sampleid fullfilepath > manifest2.csv

cat manifest.csv manifest2.csv > manifest_final.csv

qiime tools import --input-path manifest_final.csv --output-path seqs.qza --type 'SampleData[SequencesWithQuality]'
