library(BiocManager)
library(Biostrings)
library(msa)
library(rentrez)
library(pwalign)
library(ape)
library(phangorn)

# reading all of the gene information (ND2, cytochrome b)

Buteo_ND2 <-readDNAStringSet("Buteo_ND2")

Buteo_ND2

#renaming the files so that it is easier to read 
names(Buteo_ND2)<- c("albonotatus","brachyurus","galapagoensis","jamaicensis","lineatus","nitidus","platypterus","regalis","swainsoni","Harris")

Buteo_ND2

ND2_align<- msa(Buteo_ND2, method = "Muscle") #runs the msa 

ND2_align


#start cytochrome b msa 

Buteo_cytob<- readDNAStringSet("Buteo_cytob")

Buteo_cytob

names(Buteo_cytob)<- c("albonotatus","brachyurus","galapagoensis","jamaicensis","lineatus","nitidus","platypterus","regalis","swainsoni","Harris")

names(Buteo_cytob)[is.na(names(Buteo_cytob))] <- "jamaicensis"

Buteo_cytob

cytob_align<- msa(Buteo_cytob, method = "Muscle")

cytob_align

# start GAP align 

Buteo_gap<- readDNAStringSet("Buteo_GAP")

Buteo_gap

names(Buteo_gap)<- c("albonotatus","brachyurus","galapagoensis","jamaicensis","lineatus","nitidus","platypterus","regalis","swainsoni","Harris")

Buteo_gap

gap_align<- msa(Buteo_gap, method = "Muscle")

gap_align

#writing phydate for each alignment in order to use in Raxml and beast 

ND2_phy<- msaConvert(ND2_align, type = "phangorn::phyDat")

write.phyDat(ND2_phy, file = "ND2_proj_alignment.phy", format = "phylip")
#
cytob_phy<- msaConvert(cytob_align, type = "phangorn::phyDat")

write.phyDat(cytob_phy, file = "cytob_proj_alignment.phy", format = "phylip")
#
gap_phy<- msaConvert(gap_align, type = "phangorn::phyDat")

write.phyDat(gap_phy, file = "gap_proj_alignment.phy", format = "phylip")
#
ND2_phy<- msaConvert(ND2_align, type = "phangorn::phyDat")
write.phyDat(ND2_phy, file = "ND2_alignment.fasta", format = "fasta")

cytob_phy<- msaConvert(cytob_align, type = "phangorn::phyDat")
write.phyDat(cytob_phy, file = "cytob_alignment.fasta", format = "fasta")

gap_phy<- msaConvert(gap_align, type = "phangorn::phyDat")
write.phyDat(gap_phy, file = "gap_alignment.fasta", format = "fasta")
S

