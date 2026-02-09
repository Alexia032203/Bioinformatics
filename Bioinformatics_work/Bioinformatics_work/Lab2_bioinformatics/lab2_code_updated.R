pacman::p_load(tidyverse,Biostrings,msa,genepop,janitor) #installing basic packages for R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
install.packages("seqinr", dependencies = TRUE)
library(seqinr)  
BiocManager::install("pwalign")
library(pwalign)
#added all the packages used in Dr. Goren's edna example 
refs <- read.fasta("sequence.fasta", as.string = TRUE)
#based on edna example used this to read the data 

str(refs)

length(refs)

counts <- table(strsplit(refs[[1]], "")[[1]])
print(counts)

total <- sum(counts)
print(total)

percent <- (counts / total) * 100
print(percent)

print(percent["a"])
print(percent["t"])
print(percent["g"])
print(percent["c"])