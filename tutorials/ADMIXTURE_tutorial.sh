conda activate your_environment
conda install pcangsd #this might have conflicts
idev -A IBN21018 -t 02:00:00
pcangsd -beagle YOURNAME.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12
