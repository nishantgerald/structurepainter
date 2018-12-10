# source("/Users/gjohnston9/Documents/popgen/bin/RScripts/Chr_kmer_InformativeSNPs_Haploid_HMM.R") ### lab
source("/home/greg/School/popgen/bin/RScripts/Chr_kmer_InformativeSNPs_Haploid_HMM.R") ### laptop

args <- commandArgs(trailingOnly=TRUE)

# ### chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt
# homologous_filename <- args[1]

# ### all.chr.ASW_CEU_YRI.SNPs.homologous.25k.2.Q
# all_admixture_filename <- args[2]

# ### chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.Q
# admixture_filename <- args[3]

# ### ASW_Sample_IDs_Haploid.txt
# ASW_ids_filename <- args[4]

# ### chr22.Test.LA.ASW_Only.txt
# out_filename <- args[5]

# ### 5
# kmer <- as.numeric(args[6])

homologous_filename = "/home/greg/School/popgen/data/chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt"
all_admixture_filename = "/home/greg/School/popgen/data/admixture/all.chr.ASW_CEU_YRI.SNPs.homologous.25k.2.Q"
admixture_filename = "/home/greg/School/popgen/data/admixture/chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.Q"
ASW_ids_filename = "/home/greg/School/popgen/SampleIDs/ASW_Sample_IDs_Haploid.txt"
out_filename = "/home/greg/School/popgen/data/out.txt"
kmer = 5

#### Read in Strands + Select only the ASW samples ####
print("Loading Data")
ASW_Chr_Strands <- read.table(homologous_filename, head=T, sep="\t")

# eliminate duplicate rows
ASW_Chr_Strands <- subset(ASW_Chr_Strands,!(duplicated(ASW_Chr_Strands$ID)|duplicated(ASW_Chr_Strands$ID,fromLast = TRUE)))

###### Function to Create Test Chromosomes ######
Create_Test_Chromosomes <- function(Chr_Strands,Pops,Ref_IDs,Num_Switch){
  ### sample from SNPs (specifically, sample positions of SNPs)
  Test_Pos_Chromosome <- sample(Chr_Strands$POS, Num_Switch, replace=FALSE)
  # Test_Pos_Chromosome <- Chr_Strands$POS[1:Num_Switch] ### for now, pick first n rows so results can be compared with Python results
  Test_Pos_Chromosome <- sort(Test_Pos_Chromosome)

  ### get index in Chr_Strands that each SNP came from
  Test_Pos_Chr_index_stop <- sapply(1:Num_Switch, function(x){which(Chr_Strands$POS == Test_Pos_Chromosome[x])}) ### get 
  Test_Pos_Chr_index_start <- c(0,Test_Pos_Chr_index_stop)
  Test_Pos_Chr_index_stop <- c(Test_Pos_Chr_index_stop,nrow(Chr_Strands))

  Test_Chr <- numeric(nrow(Chr_Strands))
  True_Chr <- character(nrow(Chr_Strands))
  for (j in 1:length(Pops)){
    for (i in seq(j,length(Test_Pos_Chr_index_start),length(Pops))){
      Test_Chr[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]] <- Chr_Strands[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i],Ref_IDs[j]]
      True_Chr[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]] <- rep(Pops[j],length((Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]))
    }
  }
  return(list(Test_Chr,True_Chr))
}

#### Test Strand 1 ####
print("Creating Test Strands")
Known_Sample1 <- Create_Test_Chromosomes(ASW_Chr_Strands,Pops = c("Pop2","Pop1"),Ref_IDs = c("NA06984_1","NA18486_1"),Num_Switch = 4)
Test_1 <- Known_Sample1[[1]]
True_1 <- Known_Sample1[[2]]

#### Test Strand 2 ####
Known_Sample2 <- Create_Test_Chromosomes(ASW_Chr_Strands,Pops = c("Pop2","Pop1"),Ref_IDs = c("NA06984_2","NA18486_2"),Num_Switch = 6)
Test_2 <- Known_Sample2[[1]]
True_2 <- Known_Sample2[[2]]

True_LA <- data.frame(POS_Start=ASW_Chr_Strands$POS,
                      POS_End=ASW_Chr_Strands$POS,True_1,True_2)

#### Read in Admixture results for all chromosomes ####
All_ADMIX2<-read.table(all_admixture_filename,head=F)
num_pops <- ncol(All_ADMIX2)
All_ADMIX2_ordered <- All_ADMIX2[,order(colMeans(All_ADMIX2),decreasing=TRUE)]
colnames(All_ADMIX2_ordered)<-paste0("Pop",1:ncol(All_ADMIX2_ordered))

#### Read in Admixture (k=2) results for single chromosome ####
print("Loading Admixture")
ADMIX2_unordered <-read.table(admixture_filename,head=F)
Admix_chr_diff_all <- apply(ADMIX2_unordered, 2, function(x) {sum(All_ADMIX2_ordered$Pop1-x)})
ADMIX2 <- ADMIX2_unordered[,order(Admix_chr_diff_all)]
colnames(ADMIX2)<-paste0("Pop",1:ncol(ADMIX2))
rownames(ADMIX2)<-colnames(ASW_Chr_Strands)[10:ncol(ASW_Chr_Strands)]
Test_ADMIX <- rbind((table(True_LA$True_1)/nrow(True_LA)),(table(True_LA$True_2)/nrow(True_LA)))
Test_IDs <- c("Test_1","Test_2")
row.names(Test_ADMIX) <- Test_IDs
ADMIX2_with_Test <- rbind(ADMIX2,Test_ADMIX)

#### Subset ASW from full data ####
print("Subsetting - ASW Only")
ASW_ids <- read.table(ASW_ids_filename,head=F)
ASW_Only_Chr_Strands <- ASW_Chr_Strands[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as.character(ASW_ids$V1))]
ADMIX2_with_Test_ASW_Only <- ADMIX2_with_Test[c(as.character(ASW_ids$V1),"Test_1","Test_2"),]
Prob_Pops <- (colSums(ADMIX2_with_Test_ASW_Only)/nrow(ADMIX2_with_Test_ASW_Only)) #Prob(Pop1) & Prob(Pop2), respectively


#### Update ASW_Chr_Strands with Test Individuals ####
ASW_Chr_Strands_with_Test <- data.frame(ASW_Only_Chr_Strands,Test_1 = Test_1, Test_2 = Test_2)

#### Select the most informative SNPs ####
Chr_Strands_with_Test_top <- SelectInformativeSNPs(ASW_Chr_Strands_with_Test, Prob_Pops, diff_quantile = 0.1)

#### Create new data frame with k consecutive snps together (k-mers) ####
Chr_Strand_with_Test_Substrings <- CreateKmerHaplotypes(Chr_Strands_with_Test_top, kmer)

#### Create Emission Matrices for each k-mer ####
kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"), repeats=TRUE), 1, paste, collapse="")
log2_emissionmatrices <- CreateEmissionMatrix(Chr_Strand_with_Test_Substrings, ADMIX2_with_Test_ASW_Only, kmer_haplotypes)

#### Local Ancestry Inference ####
Train_with_Test_LA <- BothDirectionsLocalAncestry_Prob(Chr_Strand_with_Test_Substrings, ADMIX2_with_Test_ASW_Only, Test_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001)
Train_with_Test_LA_Both <- Train_with_Test_LA[[1]]
True_Full <- data.frame(POS=ASW_Only_Chr_Strands$POS, True_1, True_2)
True_Substrings <- merge(Chr_Strand_with_Test_Substrings[,1:3], True_Full, by.x="POS_Start", by.y="POS")
Test_LA_Reordered <- data.frame(Train_with_Test_LA_Both[,1:3], True_1=True_Substrings$True_1, Test_2=Train_with_Test_LA_Both[,4], True_2=True_Substrings$True_2)
write.table(Test_LA_Reordered, out_filename, quote=F, sep="\t", row.names = F)
out_png_filename <- gsub("txt", "png", out_filename)
PlotLocalAncestry(LocalAncestryDataFrame = Test_LA_Reordered, kmer, imagename=out_png_filename)
