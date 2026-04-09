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
library(UniprotR)
library(protti)
library(r3dmol)
# all packages loaded 

getwd() #this lets me know the working directory in order to upload files to this script 

metazoa_gene<- readDNAStringSet("metazoa_alignment.gene.fasta") #string set is now in script  

metazoa_gene # checking sequences 
names(metazoa_gene)

?grep

homo_gene<- metazoa_gene[grep("Homo_sapiens", names(metazoa_gene))] # this is used to separate the homo sapiens gene from the rest  

homo_gene # checking it 

# code was copied from lab 6 in order to change DNA to protein  

translate_code <- function(dna_seq) {
  dna_char <- as.character(dna_seq)
  dna_char <- gsub("-", "", dna_char) #removes the gaps for this to function            
  dna_char <- gsub("[^ATGCatgc]", "", dna_char) #make sure that all letters are correct 
  dna_obj <- Biostrings::DNAString(dna_char)
  Biostrings::translate(dna_obj)
}

homo_protein <- translate_code(homo_gene[[1]])
homo_protein

#following same steps as lab 10 (proteins)

protein_seq <- as.character(homo_protein) # wrote out so i can copy and paste in blastp 

protein_seq

# the accession number for the closest match to the sequence 

accession <- "P54098"

GOinfo <- GetProteinGOInfo(accession)

GOinfo #full information on three sub-ontologies 

PlotGoInfo(GOinfo) # plot has been saved as PDF















