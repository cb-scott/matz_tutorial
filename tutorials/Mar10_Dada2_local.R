#Only need to run these on the first time 
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")

library("dada2")

path <- "~/Downloads/pairedfastq/" #be sure to unzip the file 
