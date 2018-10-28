library(ChemmineOB)
library(readr)
library(Biostrings)
library(protr)
library(ChemmineR)
x = readFASTA("uniprot-yourlist_Train_target_10_24_18.fasta")
#Amino acid composition
K<-names(x)
K<-data.frame(K)
K$K <- substr(K$K, 4, 9)
names(K)<-"target_id"
COMP = t(sapply(x, extractAAC))
rownames(COMP)<-NULL
COMP<-data.frame(cbind(K, COMP))
APse = t(sapply(x, extractAPAAC))
rownames(APse)<-NULL
APse<-data.frame(cbind(K, APse))
#Depetide
DC = t(sapply(x, extractDC))
rownames(DC)<-NULL
DC<-data.frame(cbind(K, DC))
#AAindex with pc and Lag
AAA<-t(sapply(x, index =c(2:4,9,21,32,65,68,97,137,253:255,319,398:402,292), pc=5, lag=7, scale = TRUE, silent = TRUE, extractProtFP))
ALL<-data.frame(cbind(COMP,APse[,-1],DC[,-1],AAA[,-1]))
write_csv(ALL, "Train_Target_Convert_10_24_18.csv")
#Read compound data
library(readr)
Train_Final_SMILE_10_24_comp <- read_csv("Train_Final_SMILE_10_24_comp.csv")
#Merge target with compound data by target_id
Final_Train<-merge(Train_Final_SMILE_10_24_comp, ALL,by="target_id", sort=FALSE)
write_csv(Final_Train, "Training_Data_10_24_18.csv")
