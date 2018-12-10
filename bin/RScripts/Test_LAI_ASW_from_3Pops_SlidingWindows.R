source("../bin/Chr_kmer_InformativeSNPs_Haploid_HMM.R")

args <- commandArgs(trailingOnly=TRUE)
homologous_filename <- args[1]
all_admixture_filename <- args[2]
admixture_filename <- args[3]
ASW_ids_filename <- args[4]
out_filename <- args[5]
kmer <- as.numeric(args[6])

#homologous_filename <- "chr1.phase3.ASW_CEU_YRI.SNPs.homologous.txt"
#admixture_filename <- "chr1.ASW_CEU_YRI.SNPs.homologous.25k.2.Q"
#all_admixture_filename <- "../AllChr/all.chr.ASW_CEU_YRI.SNPs.homologous.25k.2.Q"
#ASW_ids_filename <- "../bin/ASW_Sample_IDs_Haploid.txt"
#out_filename <- "Testing.chr1.Test.LAI.ASW_Only.txt"
#kmer <- 5

#### Read in Strands + Select only the ASW samples ####
print("Loading Data")
ASW_Chr_Strands <- read.table(homologous_filename, head=T, sep="\t")
ASW_Chr_Strands <- subset(ASW_Chr_Strands,!(duplicated(ASW_Chr_Strands$ID)|duplicated(ASW_Chr_Strands$ID,fromLast = TRUE)))

###### Function to Create Test Chromosomes ######
Create_Test_Chromosomes <- function(Chr_Strands,Pops,Ref_IDs,Num_Switch){
  Test_Pos_Chromosome <- sample(Chr_Strands$POS, Num_Switch, replace=FALSE)
  Test_Pos_Chromosome <- sort(Test_Pos_Chromosome)
  Test_Pos_Chr_index_stop <- sapply(1:Num_Switch, function(x){which(Chr_Strands$POS == Test_Pos_Chromosome[x])})
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
#PlotLocalAncestry(LocalAncestryDataFrame = True_LA, kmer, imagename="chr1.True.Test.png")

#### Subset ASW from full data ####
print("Subseting - ASW Only")
ASW_ids <- read.table(ASW_ids_filename,head=F)
ASW_Only_Chr_Strands <- ASW_Chr_Strands[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as.character(ASW_ids$V1))]
ADMIX2_with_Test_ASW_Only <- ADMIX2_with_Test[c(as.character(ASW_ids$V1),"Test_1","Test_2"),]
Prob_Pops <- (colSums(ADMIX2_with_Test_ASW_Only)/nrow(ADMIX2_with_Test_ASW_Only)) #Prob(Pop1) & Prob(Pop2), respectively

#### Update ASW_Chr_Strands with Test Individuals ####
ASW_Chr_Strands_with_Test <- data.frame(ASW_Only_Chr_Strands,Test_1 = Test_1, Test_2 = Test_2)

#### Select the most informative SNPs ####
Chr_Strands_with_Test_top <- SelectInformativeSNPs(ASW_Chr_Strands_with_Test, Prob_Pops, diff_quantile = 0.1)

#### Local Ancestry Estimation along k windows ####
BothDirectionsListDataFrame <- replicate(kmer,data.frame())

for (window in 1:kmer){
  print(paste0("Window ",window))
  #### Create new data frame with k consecutive snps together (k-mers)####
  if (window != 1){
    Chr_Strands_top_window <- Chr_Strands_with_Test_top[-(1:(window-1)),]
  } else {Chr_Strands_top_window <- Chr_Strands_with_Test_top}
  Chr_Strand_Substrings <- CreateKmerHaplotypes(Chr_Strands_top_window,kmer)
  #### Create Emission Matrices for each k-mer ####
  kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
  log2_emissionmatrices <- CreateEmissionMatrix(Chr_Strand_Substrings,ADMIX2_with_Test_ASW_Only,kmer_haplotypes)
  #### Local Ancestry Inference for each window ####
  BothDirectionsListDataFrame[[window]] <- BothDirectionsLocalAncestry(Chr_Strand_Substrings, ADMIX2_with_Test_ASW_Only, Test_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001)
}

#### Summarize Local Ancestry Estimations across windows ####
SummaryLocalAncestryWindows <- Chr_Strands_with_Test_top[,c("POS","POS",Test_IDs)]
colnames(SummaryLocalAncestryWindows)[1:2] <- c("POS_Start","POS_End")
SummaryLocalAncestryWindows[,3:ncol(SummaryLocalAncestryWindows)] <- rep(NA, (ncol(SummaryLocalAncestryWindows)-2)*nrow(SummaryLocalAncestryWindows))

for (id in Test_IDs){
  Pop1 <- rep(0,nrow(SummaryLocalAncestryWindows))
  Pop2 <- rep(0,nrow(SummaryLocalAncestryWindows))
  for (i in 1:nrow(SummaryLocalAncestryWindows)){
    for (window in 1:kmer){
      row_id <- findInterval(SummaryLocalAncestryWindows[i,1],matrix(unlist(BothDirectionsListDataFrame[[window]][,1]),length(BothDirectionsListDataFrame[[window]][,1]),byrow=T))
      if (row_id != 0){
        if (BothDirectionsListDataFrame[[window]][row_id,id] == "Pop1"){
          Pop1[i] <- Pop1[i]+1
        } else if (BothDirectionsListDataFrame[[window]][row_id,id] == "Pop2"){
          Pop2[i] <- Pop2[i]+1
        }
      }
    }
  }
  SummaryLocalAncestryWindows[,id] <- ifelse(Pop1 == Pop2, "Unknown", ifelse(Pop1 > Pop2, "Pop1","Pop2"))
  print(paste0("Finished ", id))
}

True_Substrings <- merge(SummaryLocalAncestryWindows[,1:2],True_LA,by=c("POS_Start","POS_End"))
True_Substrings <- True_Substrings[order(True_Substrings$POS_Start),]
Test_LA_Reordered <- data.frame(SummaryLocalAncestryWindows[,1:3],True_1=True_Substrings$True_1,Test_2=SummaryLocalAncestryWindows[,4],True_2=True_Substrings$True_2)
write.table(Test_LA_Reordered, out_filename, quote=F, sep="\t", row.names = F)

#### Plot Results ####
out_png_filename <- gsub("txt", "png", out_filename)
PlotLocalAncestry(LocalAncestryDataFrame = Test_LA_Reordered, kmer, imagename=out_png_filename)
