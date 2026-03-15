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
install.packages("UniprotR")
install.packages("protti")
install.packages("r3dmol")
BiocManager::install("GenomicAlignments")
library(UniprotR)
library(protti)
library(r3dmol)
# all packages installed and loaded 

#Accession number files have been downloaded and loaded into this script
#I need to just focus on the first 5 matches 

read.delim("accs.txt")

accession_numb <- read.csv("accs.txt", header = F) 
?read.csv

accession_numb$V1[1:5] #this will make sure it is just the top 5 

accession_5<- accession_numb$V1[1:5]

?gsub

accession_5 <- gsub("SP:", "", accession_5) # this code was used to delete the SP: from the accession because that format caused the below line of code to not work 

accession_5

GetProteinGOInfo(accession_5)

Getinfo<- GetProteinGOInfo(accession_5) #assigning this to a variable 

Getinfo

PlotGoInfo(Getinfo)# plotting the results 

#handy visualization from github 

PlotGOAll(GOObj = Getinfo, Top = 10, directorypath = getwd(), width = 8, height = 5)

getwd() #checking that this went to the correct working directory  

#Q8-What are some interesting GO terms for your gene?
# the three interesting GO terms for t his gene were: 
#1. NADH dehydrogenase (ubiquinone) activity, 2. mitochondrial inner membrane, and 3. mitochondrial electron transport NADH to ubiquinone 

#Q9- this gene is the ND2 gene 

gene<- "O21398"

GetPathology_Biotech(gene)

Pathology_object <- GetPathology_Biotech(gene)

?Get.diseases

Get.diseases(Pathology_object = Pathology_object, directorypath = NULL)
# so this does not work probably because the gene not having any pathology associated with the avian version of the gene

# going to try all of this with the human version of the gene 

gene_human<- "P03891"

GetPathology_Biotech(gene_human)
# one of the diseases associated with Alzheimer disease 

Pathology_object2 <- GetPathology_Biotech(gene_human)

disease_data <- Get.diseases(Pathology_object = Pathology_object2, directorypath = NULL)

disease_data # reading the diseases associated 

# getting structural information 

fetch_uniprot(accession_5)

pdb_ids <- c("1ZMR", "2HWG") #had to use these because other ID's had nothing to retrieve  

fetch_pdb(pdb_ids)


fetch_alphafold_prediction(accession_5) # this was able to give me the information of the structures 


r3dmol(fetch_alphafold_prediction(accession_5))
# tried this but it would not load 

# had to also use the human version of the ND2 in order to get a structure in the Alpha fold 









