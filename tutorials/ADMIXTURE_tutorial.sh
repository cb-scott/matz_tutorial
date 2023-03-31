conda activate your_environment
conda install pcangsd #this might have conflicts
pcangsd -beagle YOURNAME.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12
