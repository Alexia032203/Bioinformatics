library(BiocManager)
library(msa)
library(rentrez)
library(pwalign)
library(dplyr)
library(readr)
library(seqinr)
library(Biostrings)
library(stringr)
install.packages("phangorn")
library(phangorn)

f1 <-readDNAStringSet("lineatus_ND2_1.fasta")  
f2 <-readDNAStringSet("lineatus_ND2_2.fasta")
f3 <-readDNAStringSet("lineatus_ND2_3.fasta")
f4 <-readDNAStringSet("elegans_ND2_1.fasta")  
f5 <-readDNAStringSet("elegans_ND2_2.fasta") 

merged <- c(f1, f2, f3, f4, f5)

merged #read all of them together 
# name is a little long so will change the name of each sequence 

names(merged) <- c("lineatus1", "lineatus2", "lineatus3", "elegans1", "elegans2")

merged

muscle_align <- msa(merged, method = "Muscle")

muscle_align #read it 

# another method of this function would be ClustalW

align_v2 <-msa(merged, method = "ClustalW") 

align_v2 #read it
#looks pretty similar 

#based on research it seems that Muscle is more accurate and faster (most recent tool)  

print(muscle_align)

?as.matrix
# this was advised in order to see the gaps if there are any 

align_matrix <- as.matrix(muscle_align)

dim(align_matrix)

align_matrix[, 1:80] #gives me an idea of each base comparing all sequences  


sum(align_matrix == "-")
# there is no gaps in these alignments 

alignment_length <- ncol(align_matrix)
alignment_length
#1041

all_bases <- as.vector(align_matrix)

gc_count <- sum(all_bases %in% c("G", "C"))

total_bases <- length(all_bases)

# Overall GC percentage
gc_percent <- (gc_count / total_bases) * 100

gc_percent
#47.07

#now have to convert the format to seqinr

muscle_align_seqinr <- msaConvert(muscle_align, type = "seqinr::alignment")

#calculating distance matrix 
dist_matrix <- dist.alignment(muscle_align_seqinr, "identity")

as.matrix(dist_matrix)
# these number do make sense however, I am intrested that lineatus 2 has the same distance value as elegans samples when comapred to lineatus 1 and 3   

#translating one of the sequences to amino acids (lineatus2)

translate_code <- function(dna_seq) {
  dna_char <- as.character(dna_seq)
  dna_obj <- Biostrings::DNAString(dna_char)
  Biostrings::translate(dna_obj)
 }

lineatus2_aa <- translate_code(merged[["lineatus2"]])

lineatus2_aa
# learned a lesson here you have to load seqinr after biostrings for this function to work ... this was tough 

Alignment_phydat<- msaConvert(muscle_align, type="phangorn::phyDat")


write.phyDat(Alignment_phydat, file = "alignment.fasta", format = "fasta")
