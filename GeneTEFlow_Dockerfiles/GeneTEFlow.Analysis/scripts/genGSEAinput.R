#!/home/zhongw2/afsproject/download/R/R-3.2.0/bin/Rscript
#
# this script generates input file for GSEA based on preranked gene list 
#	ranking method: signed log10 pvalue
# Author: Wenyan Zhong
# Last modified date: 07/09/2018
#

# Input files
args = commandArgs(TRUE)

DEfile <- args[1]           # deseq2 output file
samplecompare.file <- args[2]   # sample comparison file



#DEfile<-c("all.deseq2.out.txt")
#samplecompare.file<-c("../sampleinfo/samplecompare.txt")



DE<-read.table(DEfile, header=TRUE,row.names=1, check.names=F, sep="\t")
# get sample compare info
samplecompare<-read.table(samplecompare.file, header=T, sep="\t", check.names=F, stringsAsFactor=F)

## extract data for gsea only does GSEA for comparison with "Stats" of "y" in samplecomapre file 
samplecompare<-samplecompare[samplecompare$Stats %in% "y",]

totalcompare<-nrow(samplecompare)
compare<-paste(samplecompare$Group2, samplecompare$Group1, sep=".vs.")

# output samplecompare file for submitting GSEA 
output.compare.file<-c("samplecompare.gsea.txt")
compare.gsea<-paste(samplecompare$Group2, samplecompare$Group1, sep=":")
write(compare.gsea, file=output.compare.file)

for ( i in 1:totalcompare )
{
	print(compare[i])
	columns.names<-c("Log2FC", "pval", "padj")
	columns.names<-paste(compare[i], columns.names, sep=".")
    logFC.col.name<-columns.names[1]
	pval.col.name<-columns.names[2]
	fdr.col.name<-columns.names[3]
	print(columns.names)
	gsea.data<-data.frame("Name"=toupper(row.names(DE)))
	gsea.data$fcSign=sign(DE[,c(logFC.col.name)] )
#	rank by pvalue
	gsea.data$logP=-log10(DE[,c(pval.col.name)])
	gsea.data$signedlog10pval=gsea.data$logP / gsea.data$fcSign
	gsea.data<-gsea.data[,c("Name", "signedlog10pval")]
	gsea.data<-gsea.data[!is.na(gsea.data$signedlog10pval),]
	gsea.data<-gsea.data[is.finite(gsea.data$signedlog10pval),]
## rank by fdr
#	gsea.data$logQ=-log10(DE[,fdr.col.name])
#	gsea.data$signedlog10fdr=gsea.data$logQ / gsea.data$fcSign
#	gsea.data<-gsea.data[,c("Name", "signedlog10qval")]
#	gsea.data<-gsea.data[!is.na(gsea.data$signedlog10qval),]
#	gsea.data<-gsea.data[is.finite(gsea.data$signedlog10qval),]
	names(gsea.data)[1]<-"#Name"
	outputfile<-paste(compare[i], ".rnk", sep="")
	write.table(file=outputfile, gsea.data, quote=F,sep="\t",row.names=F)
}
