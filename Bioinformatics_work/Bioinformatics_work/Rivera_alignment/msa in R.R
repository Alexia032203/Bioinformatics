#Hands-on workshop for FGCU, computational biologist position
#Sequence Alignment Demo
#Date: january 26, 2026

##get rid of packages for fresh installation
remove.packages(c("Biostrings", "msa", "rentrez"))

######################
## Step 1: Set your working diretory
######################
setwd("~/Desktop/R course")

######################
## Step 2: Install the packages we will use this sessions
######################
install.packages("BiocManager")
BiocManager::install("Biostrings", force = T)
BiocManager::install("msa", force = T)
BiocManager::install("rentrez", force = T)
install.packages("rentrez")
# only way to download if you have a updated version of R


######################
## Issue: sometimes the R version and the package version don't align so you have to force install packages
######################

# !!!! OH NO !!!! This doesn't seem to work for us. Let's troubleshoot
# This code checks whehter the Bioconductor installer exists, installs it if needed, and the uses it to install
# bioinformtics packages

if (!requireNamespace("BiocManager", quietly = TRUE)) #in English: if BiocManager is not installed
+  install.packages("BiocManager", force=TRUE) #in English: Then install BiocManager (this only runs if previous line is TRUE)

BiocManager::install(c("Biostrings", "msa")) #in English: Uses BiocManager to install Biostrings and msa

#####################
# load libraries
####################
library(Biostrings)
library(msa)
library(rentrez)


#######################
## Step 3: Download DNA sequences from species we want
#######################

# Identify species you want, here we want a human, chimp, mouse, wolf, and cow
acc <- c(
  "NC_012920.1",  # Homo sapiens mitochondrion
  "NC_001643.1",  # Pan troglodytes mitochondrion
  "NC_005089.1",  # Mus musculus mitochondrion
  "NC_002008.4",  # Canis lupus familiaris mitochondrion
  "NC_006853.1"   # Bos taurus mitochondrion
         )

#get the FASTA file
fasta_text <- entrez_fetch(db = "nuccore", id = acc, rettype = "fasta", retmode = "text")
#db=databse to use, id=species we want, rettype=format of data, retmode=format to receive data

str(fasta_text) #structure of object
fasta_text

# Create a file for the sequences and saves it to your working directory
writeLines(fasta_text, "class_sequences.fasta")
#object, the file we want it in 
list.files()# in order to see what files are in teh directory  

# Read in the sequence file
dna <- readDNAStringSet("class_sequences.fasta")
dna

#change name so it's easier to read (we don't need the accession numbers anymore)
names(dna) <- c("Homo sapien", "Pan troglodytes", "Mus musculus", "Canis lupus", "Bos taurus")
dna

#######################
## Step 4: Clean your data, this code replaces ambiguities with N's (languages don't like \n, they mean something specific)
#######################
dna.clean <- replaceAmbiguities(dna, new="N")

dna.clean

#######################
## Step 5: Align your cleaned sequences using MUSCLE, Progressive method
# MUSCLE builds an alignment by progressively aligning sequences, then repeatedly refining it to improve accuracy, this is the default method in most pipeline
#######################
aln <- msa(dna.clean, method = "Muscle")
class(aln)
aln







## If there is extra time, we can make a phylogeny ##

install.packages("ape")




#######################
## Step 6: Make a tree
#######################
aligned <- as.DNAbin(aln)  #makes the seuqences useful for downstream analyses
aligned

d <- dist.dna(aligned, model = "K80") 
#convers your alignment into a table of "how evolutionary different is each species from every other species"
#method using called Kimura 2-paramter (transitions and transversions occur at different rates)

#######################
# Food for thought
#######################
#Building a tree based on pure similarity, can you figure out what the tree will look like just by using the distance matrix values?
d


hc <- hclust(as.dist(d), method = "average")   # UPGMA-style clustering, UPGMA build a tree by repeatedly grouping together the most similar sequences, step-by-step.
tree <- as.phylo(hc) #forces the heuristic clustering into a phylogeny

class(tree) #if the previous code worked, this should say "Phylo" or "MultiPhylo" if many trees
tree


#######################
## Step 7: Plot a tree
#######################
plot(tree, main = "My First Phylongey")
















