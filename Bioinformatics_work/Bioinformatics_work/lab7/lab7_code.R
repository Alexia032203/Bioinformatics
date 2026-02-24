library(BiocManager)
library(msa)
library(rentrez)
library(pwalign)
library(dplyr)
library(readr)
library(seqinr)
library(Biostrings)
library(stringr)
library(phangorn)

ref_61800 <- readDNAStringSet("Contopus_sordidulus_61800_ref.fasta") # this is Contopus sordidulus reference 

ref_61800


ref_768<- readDNAStringSet("Contopus_virens_768_ref.fasta") # this is Contopus virens reference  

ref_768

unknown_files<- c("Species_identifier_consensus_curated_90763.fa_1.fa",
                  "Species_identifier_consensus_curated_91911.fa.fa",
                  "Species_identifier_consensus_curated_91929.fa")
unknown_files

unknown_seqs <- readDNAStringSet(unknown_files)

names(unknown_seqs) <- c("Unknown_90763", "Unknown_91911", "Unknown_91929") # rename for easier comprehension 

unknown_seqs

#doing pairwise alignment based on Dr. Schulte code on this similar topic 

combined_seqs <- c(ref_61800, ref_768, unknown_seqs)

combined_seqs
names(combined_seqs)

aln_90763_1<- pairwiseAlignment(combined_seqs[3], combined_seqs[1])
pid(aln_90763_1) #93.97

aln_90763_2<- pairwiseAlignment(combined_seqs[3], combined_seqs[2])
pid(aln_90763_2) #87.17 

# 90763 is a pretty strong match for Contopus sordidulus 

aln_91911_1<- pairwiseAlignment(combined_seqs[4], combined_seqs[1])
pid(aln_91911_1) #91.05

aln_91911_2<- pairwiseAlignment(combined_seqs[4], combined_seqs[2])
pid(aln_91911_2) # 99.09 

# 91911 is avery strong match to Contopus virens 

aln_91929_1<- pairwiseAlignment(combined_seqs[5], combined_seqs[1])
pid(aln_91929_1) #96.94

aln_91929_2<- pairwiseAlignment(combined_seqs[5], combined_seqs[2])
pid(aln_91929_2)#86.25 

# 91929 is a pretty strong match for Contopus sordidulus 






