install.packages("BiocManager")
BiocManager::install("Biostrings", force = T)
BiocManager::install("msa", force = T)
BiocManager::install("rentrez", force = T)
install.packages("rentrez")
setwd("~/Bioinformatics_work/Schulte_alignment")
library(BiocManager)
library(Biostrings)
library(msa)
library(rentrez)
library(pwalign)

seqs<- readDNAStringSet("Schulte_pairwiseAlignmentSequences.fasta")

seqs
# read the sequences first 

names(seqs)
# large mouth bass, small mouth bass, 2 unknown
#COI used to identify the species 

aln<- pairwiseAlignment(seqs[1], seqs[2])
aln
str(aln)

writePairwiseAlignments(aln)

nmatch(aln)
# 1432

nmatch(aln)/width(seqs[1])* 100

pid(aln)

aln2<- pairwiseAlignment(seqs[3], seqs[1])
#match percentage 92.33 (distance) 

pid(aln2)
# total match percentage between a large mouth and unknown 96.78 

aln3<- pairwiseAlignment(seqs[3], seqs[2])

pid(aln3)
# total match percentage between a large mouth and unknown 92.26

# the unknown species is closer to a large mouth bass 
# this is a different species but recently so (this is a Florida bass)
