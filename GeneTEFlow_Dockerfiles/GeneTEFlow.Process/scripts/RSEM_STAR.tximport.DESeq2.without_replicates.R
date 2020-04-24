library(readr)
library(plyr)

args = commandArgs(TRUE)

sampleinfoFile<-args[1]   #sample information file
samplecompareFile<-args[2] # group comparison design
compare.group<-args[3]  # comparison group label, an RNAseq may have subgroups for comparison
tx2genefile<-args[4] #the file for tximport


sampleinfo<-read.table(sampleinfoFile, header=TRUE, sep="\t", stringsAsFactors=F )
samplecompare<-read.table(samplecompareFile, header=T, sep="\t", check.names=F, stringsAsFactors=F)


samplecompare=na.omit(samplecompare)

comparegroup<-c(samplecompare$Group2, samplecompare$Group1)
comparesampleinfo<-sampleinfo[sampleinfo$SampleGroup %in% comparegroup,]

totalcomp<-nrow(samplecompare)

# subset data
files <- file.path(paste (comparesampleinfo$FastqID, "RSEM_Output.genes.results",sep = ".", collapse = NULL))
names(files) <- comparesampleinfo$FastqID
all(file.exists(files))




tx2gene <- read_csv(tx2genefile, col_names = FALSE)
colnames(tx2gene)=c("TXNAME ","GENEID")




library(tximport)
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


#txi has three attributes:
     #1) txi$abundance: which from "TPM" column of *.RSEM_Output.genes.results
     #2) txi$counts: which from "expected_count" column of *.RSEM_Output.genes.results
     #3) txi$length: which from "effective_length" column of *.RSEM_Output.genes.results




# filter data
#rs=rowSums(data)
#use =  (rs > 10)
print("Before filtering:")
print(summary(txi))
rowMax<-apply( txi$counts, 1, max)
use =  (rowMax >=10 )
txi$counts=txi$counts[use,]
txi$abundance=txi$abundance[use,]
txi$length=txi$length[use,]
print("After filtering:")
print(summary(txi))




# if length ==0 in some cases, Return an error - all(lengths > 0) is not TRUE.
# first define a function that counts the number of zero elements in a vector, then test the matrix columns
count0 <- function(x){length(which(x==0))}
rsem.zeros <- apply(txi$length,2,count0)
print("The number of zero elements: ")
print(rsem.zeros)

#fix the error for DESeq2
txi$length[txi$length == 0] <- 1

library(DESeq2)


fileTable=as.data.frame(row.names=comparesampleinfo$FastqID,cbind(files, comparesampleinfo$FastqID))
colnames(fileTable)=c("Filename","FastqID")


colData<-data.frame(row.names=comparesampleinfo$FastqID, FastqID=comparesampleinfo$FastqID)
colData<-join(colData, fileTable, by="FastqID")
row.names(colData)<-colData$FastqID

colData<-join(colData, comparesampleinfo, by="FastqID")
colData$SampleGroup=factor(colData$SampleGroup)
row.names(colData)<-colData$FastqID




dds<-DESeqDataSetFromTximport(txi,colData = colData,design = ~ SampleGroup)



dds <- DESeq(dds)



norm.counts<-counts(dds, normalized=TRUE)
samplename<-comparesampleinfo$SampleName
sampleid<-comparesampleinfo$FastqID
count.m.id<-colnames(norm.counts)
count.m.name<-samplename[match(count.m.id,sampleid)]
colnames(norm.counts) <- count.m.name

tmpres<-results(dds)
tmpres<-as.data.frame(tmpres)
res.all<-data.frame(row.names=row.names(tmpres))
for ( i in 1:totalcomp )
{
        group1<-samplecompare$Group1[i]
        group2<-samplecompare$Group2[i]
        groups<-c(group1, group2)
        res<-results(dds, contrast=c("SampleGroup", group2, group1))
        res<-as.data.frame(res)

        baseMeanPerLvl <- sapply( c(group1, group2), function(lvl) ( counts(dds,normalized=TRUE)[,dds$SampleGroup == lvl] ) )
        outresult<-data.frame(baseMeanPerLvl[,c(1,2)], res[,c(2,5,6)])
        names(outresult)[1:2]<-paste(names(outresult)[1:2], ".Mean", sep="")
        name.prefix<-paste(group2, ".vs.", group1, sep="")
        names(outresult)[3:5]<-paste(name.prefix, c( "Log2FC", "pval", "padj"), sep=".")
        res.all<-cbind(res.all, outresult)
}
output<-cbind(res.all, norm.counts)
outputfile<-paste(compare.group, ".deseq2.RSEM_STAR.out.txt", sep="")
write.table(file=outputfile, cbind(GeneSymbol=row.names(output), output), row.name=F, sep="\t", quote=F)


