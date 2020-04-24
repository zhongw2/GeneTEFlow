library(xlsx)

args = commandArgs(TRUE)

sampleinfoFile <- args[1]   #sample information EXCEL file
file.extention <- args[2]   #string for clip FastqID
                            #file.extention<-c("_R\\d+.fastq.gz")
                            #Or 
                            #file.extention<-c("_\\d+.clipped.fastq.gz")

sample.manifest.sheetname <- args[3]  #sample.manifest sheetname
samplecompare.sheetname <- args[4]    #samplecompare sheetname

sample.manifest <- read.xlsx(sampleinfoFile, sheetName=sample.manifest.sheetname)
sample.manifest$FastqID<-gsub(file.extention, "", sample.manifest$FastqID, perl=T)
#sample.info<-unique(sample.manifest[, c(1,3:6)])
sample.info<-unique(sample.manifest[, 1:6])



sampleinfo.file<-c("sampleinfo.txt")
write.table(file=sampleinfo.file, sample.info, sep="\t", row.names=F, quote=F)


samplecompare.file<-c("samplecompare.txt")
sample.compare<- read.xlsx(sampleinfoFile, sheetName=samplecompare.sheetname)
write.table(file=samplecompare.file,sample.compare, sep="\t", row.names=F, quote=F)

