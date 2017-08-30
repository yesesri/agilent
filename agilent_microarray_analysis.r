#set your working path 
setwd("/home/ycheruku/agilent")
#ftp link to download R-packages
source("http://bioconductor.org/biocLite.R") 
#download GEOquery package to retrive datasets from geo database
biocLite("GEOquery") 
#package to do extract expression data from Agilent sequenced microarray data and to do DE analysis
biocLite("limma") 
#call the installed packages
library(limma) 
library(GEOquery) 
#download the completedataset of the study from GEO database.
getGEOSuppFiles("GSE36875") 
##unzip and make sure all your samples are available and remove remainaning samples####
setwd("/home/ycheruku/agilent/GSE36875/")
dat=read.maimages(files=dir(), source="agilent")
#extract expression data from all your sample files in your current directory
raw_BGcorrected =backgroundCorrect(dat, "normexp", offset=50) # normalize the expression across the samples
#Then normalize and log-transformed the data. 
raw_BGandNormalized = normalizeBetweenArrays(raw_BGcorrected,method="quantile")
#Finally calculate the average intensity values from the probes of each gene.   
raw_aver = avereps(raw_BGandNormalized,ID=raw_BGandNormalized$genes$ProbeName)
##### Differential Expression Analysis ############
setwd("/home/ycheruku/agilent")
##load package ##
library(limma)
##read phenotypefile ## ## defines the categories in the data ##
pheno =read.table("pheno_file.txt",row.names=1)
## load the expression file ##
data<- raw_aver
## prepare a design matrix ##
Group<-factor(pheno$V2,levels=levels(pheno$V2))
design<-model.matrix(~0+Group)
colnames(design)<-levels(Group)
## fit the linear model ##
fit = lmFit(data,design)  ### make sure your colnames in data and rownames in design are same and equal#####
## perform the statistical Analysis ##
fit = eBayes(fit)
## report the top 20 Differentially Expressed genes ##
x<-topTable(fit,number=20,coef =2)
## write the Differentially Expressed genes to the output file ## 
write.table(x,file="DE_result.txt",sep="\t")
