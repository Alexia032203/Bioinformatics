if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

pacman::p_load(tidyverse,Biostrings,msa,genepop,janitor,seqinr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pwalign", ask = FALSE, update = FALSE)
library(pwalign)

#added all the packages used in Dr. Goren's edna example
#Escherichia coli read 
refs_ecoli <- read.fasta("sequence.fasta", as.string = TRUE)
#based on edna example used this to read the data 

str(refs_ecoli)

length(refs_ecoli)

counts_ecoli <- table(strsplit(refs_ecoli[[1]], "")[[1]])
print(counts_ecoli)

total_ecoli <- sum(counts_ecoli)
print(total_ecoli)

percent_ecoli <- (counts_ecoli / total_ecoli) * 100
print(percent_ecoli)

print(percent_ecoli["a"])
print(percent_ecoli["t"])
print(percent_ecoli["g"])
print(percent_ecoli["c"])

gc_percent <-percent_ecoli["g"]+ percent_ecoli["c"]  
print(gc_percent)

# now I can read thermoplasma volcanium

refs_thermo <- read.fasta("thermo_vol.fasta", as.string = TRUE)
str(refs_thermo)

length(refs_thermo)

counts_thermo <- table(strsplit(refs_thermo[[1]], "")[[1]])
print(counts_thermo)
# a = 478,011, c = 315,483, g = 317,147, t = 474,163

total_thermo <- sum(counts_thermo)
print(total_thermo)
# total = 1,584,804

percent_thermo <- (counts_thermo / total_thermo) * 100
print(percent_thermo)

print(percent_thermo["a"])
# 24.62
print(percent_thermo["t"])
#24.59
print(percent_thermo["g"])
#25.37
print(percent_thermo["c"])
#25.42

gc_percent_thermo <-percent_thermo["g"]+ percent_thermo["c"]  
print(gc_percent_thermo)
#50.79