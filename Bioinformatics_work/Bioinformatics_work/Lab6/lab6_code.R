library(BiocManager)
library(Biostrings)
library(msa)
library(rentrez)
library(pwalign)
library(dplyr)
library(readr)
library(seqinr)

f1 <-readDNAStringSet("lineatus_ND2_1.fasta")  
f2 <-readDNAStringSet("lineatus_ND2_2.fasta")
f3 <-readDNAStringSet("lineatus_ND2_3.fasta")
f4 <-readDNAStringSet("elegans_ND2_1.fasta")  
f5 <-readDNAStringSet("elegans_ND2_2.fasta") 

merged <- c(f1, f2, f3, f4, f5)

merged #read all of them together 







