library(splitstackshape)

args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]

Chr <- read.table(filename,head=T,comment.char = "#")
Chr_Strands <- data.frame(cSplit(Chr, colnames(Chr)[10:ncol(Chr)], sep = "|", direction = "wide", fixed = TRUE,
                                 drop = TRUE, stripWhite = TRUE, makeEqual = NULL, type.convert = TRUE))
allele_filename <- gsub("recode.vcf", "homologous.txt", filename)
write.table(Chr_Strands,allele_filename,col.names = T,row.names=F,quote=F,sep="\t")
Chr_Strands[,10:ncol(Chr_Strands)] <- apply(Chr_Strands[,10:ncol(Chr_Strands)],2,function(x){paste(x,x,sep="|")})
homologous_filename <- gsub("recode.vcf", "homologous.vcf", filename)
write.table(Chr_Strands,homologous_filename,col.names = T,row.names=F,quote=F,sep="\t")