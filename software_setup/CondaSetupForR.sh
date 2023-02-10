#alternative tutorial from bioconda
https://bioconda.github.io/user/install.html

##HOW I GOT IT TO WORK
####################################################
########   DOWNLOAD ANACONDA #######################
####################################################

#download anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
#Run the file to install from source
#IMPORTANT: Put in WORK, NOT HOME! It is too big to house in your home directory
#I recommend creating a software subdirectory in $WORK or $STOCKYARD for it.
#If you put it in stockaryd, I think you should be able to use the environments
#on any supercomputing cluster.

#When you go to install anaconda, it will ask you what you want the library path to be
#Be sure to give it the path to your WORK/STOCKYARD folder or it will not be usable.

#I believe conda will automatically update your .bashrc
#NOTE that after you install anaconda it will supersede any of the modules you
#have loaded in TACC if there is *also* an anaconda installation.
#In other words, you can still load and run TACC modules, but if you have the
#same software downloaded through conda, that version will take preference.

#If you change your mind and no longer want to run anaconda, you will have to
#uninstall it to undo these changes. To do so:

#https://docs.anaconda.com/anaconda/install/uninstall/

#After installing, add channels you will use to install packages.
#Basically, this tells anaconda where to find the packages you will want to use
#Be sure to do this in this order!! This sets the priority of the different channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


####################################################
########  SET UP YOUR ENVIRONMENT ##################
####################################################
#The big catch after you download anaconda is that you will need
# to make your first environment. You cannot download any packages until you do this

#create a new conda environemnt to house my r packages etc
#My environment will be called BleachModel for my modeling.
conda create --name BleachModel

#activate the environment
conda activate BleachModel
#note: conda activate will also let you switch between environments later on
#This way if you wanted to run two different softwares,

#check that the environment is active,
#it should have a * next to it if it is active
conda env list



####################################################
########  INSTALL PACKAGES ########################
####################################################

#which packages are active in your environment?
#Basically check here to see if what you want is downloaded
conda list


#Install new packages to conda
#sometimes these packages have slightly different names than what you would expect.
#Just google: conda "mypackage" and it will give you the exact command to use
#If you want, you can install may packages at once. After you get the hang of this
#it is much faster to run as a job on tacc than on the login node.
conda install mypackage

#EXAMPLE, instead of just "r", you will:
conda install -c conda-forge r-base=4.1.2
#here I tell conda that I specifically want R version 4.1.2 from the channel
#conda-forge. Channels are similar to CRAN repositories in R.

#EXAMPLE, if you wanted a local bowtie2
#This one is trickier, because bowtie2 is from the "bioconda" channel.
#We don't have to specify the channel, becuase we already added it to our
#environment on line 33.
#Further, unlike the r installation, since there is only one
#channel that has bowtie2, we don't need a ""-c" argument
conda install bowtie2

####################################################
#############  RUN PACKAGES ########################
####################################################

#The best part? After installing any software, it is automatically
#added to your path by conda. Now if I want to run bowtie2, instead of
#having to module load bowtie, I can simply

bowtie2 'my arguments here'

#(Last I checked - on stampede) This also works when submitting jobs to the compute nodes -
#you should not have to load any modules in your script file before submitting.
