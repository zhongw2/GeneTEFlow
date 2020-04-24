args = commandArgs(TRUE)

sampleinfo.file <- args[1]   #sample information txt file sampleinfo.txt
tpm.file <- args[2]         #TPM file all.sample.TPM.genes.results
count.file <- args[3]         #Counts file all.sample.Counts.genes.results







#tpm.file<-c("all.sample.TPM.genes.results")
#count.file<-c("all.sample.Counts.genes.results")
#sampleinfo.file<-c("../sampleinfo/sampleinfo.txt")

sampleinfo<-read.table(sampleinfo.file, header=T, sep="\t", check.names=F, stringsAsFactor=F)
sampleid<-sampleinfo$FastqID
samplename<-sampleinfo$SampleName

tpm<-read.table(tpm.file, header=T, check.names=F, stringsAsFactors=F, row.names=1, sep="\t")
names(tpm)<-gsub(".genes.results", "", names(tpm), perl=T)

count<-read.table(count.file, header=T, check.names=F, stringsAsFactors=F,row.names=1, sep="\t")
names(count)<-gsub(".genes.results", "", names(count), perl=T)

tpm.rda.file<-c("all.sample.tpm.rda")
count.rda.file<-c("all.sample.counts.rda")

saveRDS(file=tpm.rda.file, tpm)
saveRDS(file=count.rda.file, count)

tpmid<-names(tpm)
#tpmname<-samplename[match(tpmid, sampleid)]
tpmname<-samplename[match(tpmid, samplename)]
names(tpm)<-tpmname

countid<-names(count)
#countname<-samplename[match(countid, sampleid)]
countname<-samplename[match(countid, samplename)]
names(count)<-countname

tpm.txt.file<-c("all.sample.tpm.txt")
write.table(file=tpm.txt.file, cbind(GeneSymbol=row.names(tpm), tpm), row.names=F, quote=F, sep="\t")
count.txt.file<-c("all.sample.counts.txt")
write.table(file=count.txt.file, cbind(GeneSymbol=row.names(count), count),row.names=F, quote=F, sep="\t")
