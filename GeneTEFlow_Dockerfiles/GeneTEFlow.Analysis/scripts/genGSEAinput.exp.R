#!/home/zhongw2/afsproject/download/R/R-3.2.0/bin/Rscript
#
# this script generates input file for GSEA based on tpm value 
# Author: Wenyan Zhong
# Last modified date: 07/09/2018
#

## input files
args = commandArgs(TRUE)

TPMfile <- args[1]              # TPM rData file all.sample.tpm.rda
sampleinfofile <- args[2]       # sampleinfo file
samplecomparefile <- args[3]   # sample comparison file



#TPMfile<-c("all.sample.tpm.rda")
#sampleinfofile<-c("../sampleinfo/sampleinfo.txt")
#samplecomparefile<-paste("../sampleinfo/samplecompare.txt")

sampleinfo<-read.table(sampleinfofile, header=T, sep="\t", check.names=F, stringsAsFactors=F)
samplecompare<-read.table(samplecomparefile, header=T, sep="\t", check.names=F, stringsAsFactors=F)

## extract data for gsea
samplecompare<-samplecompare[samplecompare$Stats %in% "y",]

tpm<-readRDS(TPMfile)

fastqid<-sampleinfo$FastqID
samplename<-sampleinfo$SampleName
tpmid<-names(tpm)
#tpmname<-samplename[match(tpmid, fastqid)]
#names(tpm)<-tpmname

totalcompare<-nrow(samplecompare)
compare<-paste(samplecompare$Group2, samplecompare$Group1, sep=".vs.")
# output samplecompare file for GSEA 
output.compare.file<-c("samplecompare.gsea.exp.txt")
compare.gsea<-paste(samplecompare$Group2, samplecompare$Group1, sep=":")
write(compare.gsea, file=output.compare.file)

for ( i in 1:totalcompare )
{
	group1<-samplecompare$Group1[i]
	group2<-samplecompare$Group2[i]
	group1.samplename<-sampleinfo$SampleName[sampleinfo$SampleGroup %in% group1]
	group1.numsample<-length(group1.samplename)
	group2.samplename<-sampleinfo$SampleName[sampleinfo$SampleGroup %in% group2]
	group2.numsample<-length(group2.samplename)
	data<-tpm[, c(group1.samplename, group2.samplename)]
	data.max<-apply( data, 1, max) 
#	print(names(data))
	use<-(data.max> 0)
	data.f<-data[use,]
	numgene<-nrow(data.f)
	numsample<-ncol(data.f)
	data.f<-data.frame("GeneSymbol"=row.names(data.f),"Description"=row.names(data.f), data.f )
#	print(data.f[1,])
#	print(numgene)
#	print(numsample)
	gctline1<-c("#1.2")
	gctline2<-c(numgene, numsample)
	gctfile<-paste(compare[i], "gct", sep=".")
#	print(gctfile)
	write(gctline1, gctfile, ncolumns=1 )
	write(gctline2, gctfile, ncolumns=2, sep="\t", append=T)
	write.table(data.f, gctfile, quote=FALSE, sep="\t", row.names=FALSE, append=T) 

	clsline1<-c(numsample, 2, 1)
	clsline2<-c("#", group1, group2)
	clsline3<-c(rep(group1,group1.numsample), rep(group2,group2.numsample))
	clsfile<-paste(compare[i], "cls", sep=".")
	write(clsline1, clsfile, ncolumns=3, sep="\t")
	write(clsline2, clsfile, ncolumns=3, sep="\t", append=T)
	write(clsline3, clsfile, ncolumns=numsample, sep="\t", append=T)
}

