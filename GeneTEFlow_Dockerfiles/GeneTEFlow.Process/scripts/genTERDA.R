args = commandArgs(TRUE)

sampleinfo.file <- args[1]   #sample information txt file sampleinfo.txt
count.file <- args[2]         #Counts file all.sample.Counts.genes.results







#tpm.file<-c("all.sample.TPM.genes.results")
#count.file<-c("all.sample.Counts.genes.results")
#sampleinfo.file<-c("../sampleinfo/sampleinfo.txt")

sampleinfo<-read.table(sampleinfo.file, header=T, sep="\t", check.names=F, stringsAsFactor=F)
sampleid<-sampleinfo$FastqID
samplename<-sampleinfo$SampleName


count<-read.table(count.file, header=T, check.names=F, stringsAsFactors=F,row.names=1, sep="\t")
names(count)<-gsub(".TE.results", "", names(count), perl=T)


tpm.rda.file<-c("all.sample.tpm.rda")

count.rda.file<-c("all.sample.counts.rda")

#saveRDS(file=tpm.rda.file, tpm)
saveRDS(file=count.rda.file, count)


countid<-names(count)
#countname<-samplename[match(countid, sampleid)]
countname<-samplename[match(countid, samplename)]
names(count)<-countname

############
#convert read counts to expression values via TPM and return these values
#https://support.bioconductor.org/p/91218/
###########
counttmp=as.data.frame(rownames(count))
colnames(counttmp)=c("fullstring")
counttmp$geneLength = 0



geneLengthfunction <- function(x) {
 start = unlist(strsplit(x,split='|',fixed=TRUE))[2]
 end = unlist(strsplit(x,split='|',fixed=TRUE))[3]
 length = abs(as.integer(end) - as.integer(start) ) + 1
 return(length)
}


for (row in 1:nrow(counttmp)) {
     counttmp[row,2]=geneLengthfunction(as.vector(counttmp[row,1]))
}

gene.length = counttmp[,2]

tpm.temp <- count / gene.length
tpm <- t( t(tpm.temp) * 1e6 / colSums(tpm.temp) )
names(tpm)<-countname
saveRDS(file=tpm.rda.file, tpm)


count.txt.file<-c("all.sample.counts.txt")
write.table(file=count.txt.file, cbind(TEid=row.names(count), count),row.names=F, quote=F, sep="\t")


tpm.txt.file<-c("all.sample.tpm.txt")
write.table(file=tpm.txt.file, cbind(TEid=row.names(tpm), tpm), row.names=F, quote=F, sep="\t")


