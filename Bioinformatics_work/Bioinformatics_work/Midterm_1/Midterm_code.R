library(ape)
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

# downloaded all packages that may be used for this code 

seqs<-readDNAStringSet("sequences.fasta") #read all sequences 

seqs

length(seqs) # twenty sequences total, looks good 

alignment<- msa(seqs, method = "ClustalW") # chose to do clustalW becuase the sequences aren't very long so this should work well   

alignment

# in order to check if this alignment is good or bad, I will check percent identity

# have to convert this to be able to do this step 

?as.DNAbin# this function is in ape and works best for converting msa for pid   

alignment_dna<- as.DNAbin(alignment)

pid_matrix<- dist.dna(alignment_dna, model = "raw", pairwise.deletion = T)

pid_matrix # alot of these numbers are 0.000 

percent_identity<- 1- as.matrix(pid_matrix) # converted numbers to make reading them easier 1.00 indicates that the alignment is identical 

percent_identity #lowest percentage was 0.986 which is still pretty similar so this looks like a good alignment 
#Homo_sapiens_6 is a sequence of interest


#calculating the consensus sequence

consensus_seq<- consensusString(alignment)

consensus_seq


# calculating GC content 

consenus_DNA<- DNAString(consensus_seq) 

gc_count<- letterFrequency(consenus_DNA, letters =c("G", "C"), as.prob = T)

gc_percent<- sum(gc_count)*100

gc_percent # 51.56%

#question 5
# based on the percent identity the sequence that seems to be different is Homo_sapiens_6 
#based on the align_matrix, there is a - at the first nucleotide rather than the A in all the rest of the sequences
# there are also two substitutions an A replacing a C and a G replacing A 

#question 6 
# I pasted Homo_sapiens_6 into BLAST, the gene is mutant hemoglobin beta chain gene (HBB), the accession number to the best match is AY356351.1  

#question 7, translating sequence 6 to proteins 

HBB_seq<-readDNAStringSet("seqdump.txt") # so i tried to do this with the 6th sequence from the original seq but it came up with errors and the sequence had multiple stop codons and gaps   
#I downloaded the sequence from blast that was a 100% match and read it in here for translation    

HBB_seq

translate_code <- function(dna_seq) {
  dna_char <- as.character(dna_seq)
  dna_obj <- Biostrings::DNAString(dna_char)
  Biostrings::translate(dna_obj)
}

protein_seq2<- translate_code(HBB_seq)

protein_seq2

?AAStringSet

protein_HBB<- AAStringSet(protein_seq2)
names(protein_HBB)<- "mutant_protein"

Biostrings::writeXStringSet(protein_HBB, filepath = "Homo_sapiens_6_protein.fasta")

#now this is a fasta file, could not do phydat because its just one sequence 

#question 8, used BLASTp and put this sequence in a query to find the closest match, used non-redundant protein sequences   
# the closest match accession number is XP_011830555.1

#question 9, in OMIM I searched for HBB
#this gene is associated with sickle cell diseases 
# based on the consisted matches with HBB I can assume that this individual does have this disease













