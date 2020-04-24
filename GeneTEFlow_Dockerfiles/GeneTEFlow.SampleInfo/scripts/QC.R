library(lattice)
library(ggplot2)
library(reshape2)
#library(gplots)
library("RColorBrewer")

args = commandArgs(TRUE)

sampleinfo.file <- args[1]   #sample information txt file sampleinfo.txt
data.file <- args[2]         #TPM data all.sample.tpm.rda



#data.file<-c("all.sample.tpm.rda")  # TPM data
#sampleinfo.file<-c("sampleinfo.txt")  # sample information file

#get exp values 
data.exp<-readRDS(data.file)

names(data.exp) <- gsub(".genes.results", "", names(data.exp), perl=T)
# get sample information
sampleinfo<-read.table(sampleinfo.file, header=TRUE, sep="\t",stringsAsFactors=F )

# get sample annotation
#sampleid<-sampleinfo$SampleID
#sampleid<-sampleinfo$FastqID
sampleid<-sampleinfo$SampleName
samplename<-sampleinfo$SampleName
samplegroup<-sampleinfo$SampleGroup
samplemodel<-sampleinfo$CellLine
tpmid<-names(data.exp)
tpmname<-samplename[match(tpmid, sampleid)]
tpmgroup<-samplegroup[match(tpmid,sampleid)]
tpmmodel<-samplemodel[match(tpmid, sampleid)]
## sample TPM summary
data.stat <-  data.frame(row.names=tpmid)
data.stat$count_nonzero <- apply ( data.exp, 2, function (x) length(x[x!=0] ))
data.stat$per_nonzero <- apply ( data.exp, 2, function (x) round(length(x[x!=0])/length(x)*100, digits=2) )

data.stat$count_mt100 <- apply ( data.exp, 2, function (x) length(x[x>=100] ))
data.stat$per_mt100 <- apply ( data.exp, 2, function (x) round(length(x[x>=100])/length(x) *100, digits=2 ) )

data.stat$reads_median <- apply (data.exp, 2, function (x) round(quantile(x, probs=0.5), digits=0) )
data.stat$reads_75th <- apply (data.exp, 2, function (x) round(quantile (x, probs =0.75 ),digits=0) )
data.stat$reads_90th <- apply (data.exp, 2, function (x) round(quantile (x, probs =0.90 ), digits=0) )

write.table(file="all.sample.TPM.stat.txt", cbind(Sample=row.names(data.stat), data.stat), row.names=F, col.names=T, sep="\t", quote=T)

## generate QC Plots
sample.col<-rainbow(ncol(data.exp), start=0.1, end=0.9)
stacked.data.exp<-stack(data.exp)

graph.file<-c("Sample.QC.pdf")
pdf(file=graph.file, width=11, height=8)
# Box Plot
box.title<-c("Box Plot of All Samples")
par(mar=c(15, 5,5,5))
boxplot(log2(data.exp+1), border="blue", main=box.title, col=sample.col, boxwex=0.7, horizontal=FALSE, 
		ylab="log2(TPM+1)", cex.main=1.5, las=2, cex=1.0, axes=F, cex.lab=1.5 )
axis(2, las=2, cex.axis=1.2)
axis(1, at=seq(1, ncol(data.exp)), label=names(data.exp), las=2, cex.axis=1.2)
box( col='black')

# Histogram Distribution
hist.title<-c("Expression Distribution of All Samples" )

#densityplot(~values, groups=ind,  data=stacked.data.exp,  plot.points = FALSE, ref=TRUE,
densityplot(~log2(values+1), groups=ind,  data=stacked.data.exp,  plot.points = FALSE, ref=TRUE,
            xlab=list(label="log2(TPM+1)", cex=1.2),
            ylab=list(label="Density", cex=1.2),
            main=list(label=hist.title, cex=1.2),
            scale = list( x= list(cex=1.2), y=list(cex=1.2)),
            auto.key=list( col=sample.col, space="right", cex=0.8, size=0.5 ),
            par.settings = list (superpose.line=list(col=sample.col ))
)


# Sample Correlation 
cor.title<-c("Sample Correlation" )
cor.melt.data<-melt(cor(data.exp))
#cor.melt.data<-melt(cor(log2(data.exp+1)))
names(cor.melt.data)[3]<-c("r")
meas <- as.character(unique(cor.melt.data$Var2))
q<-qplot(x=Var1, y=Var2, data=cor.melt.data, 
	fill=r, geom="tile", xlab="", ylab="", main=cor.title)
q+scale_fill_gradient(low="green", high="red", limits=c(0.5, 1) ) + 
 theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, colour="black"), 
		axis.text.y=element_text(colour="black") ) + 
# geom_text(label=c(round(cor.melt.data$r,2)), colour="white", size=3) +
 scale_x_discrete(limits=meas[1:length(meas)]) +  
 scale_y_discrete(limits=meas[1:length(meas)]) 

# sample correaltion with scatter plot 
panel.cor <- function(x, y, digits=3, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r.p <- abs(cor(x, y, method='pearson'))
    txt.p <- format(c(r.p, 0.123456789), digits=digits)[1]
    txt.p <- paste("r(pearson)=", txt.p, sep="")
    if(missing(cex.cor)) cex <- 0.5/strwidth(txt.p)
    test.p <- cor.test(x,y, method='pearson')
    # borrowed from printCoefmat
    Signif.p <- symnum(test.p$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
    pval.p<-paste("p(pearson)=", format(test.p$p.value, scientific=TRUE, digits=
2), sep="")
    text(0.5, 0.6, txt.p, cex=cex)
    text(0.5, 0.5, pval.p, cex=cex)
    text(.8, .9, Signif.p, cex=cex, col=2)
}

#comp<-log2(data.exp[,c(21:25)]+1)
#for ( i in c(0,1))
#{
#	index<-i*5 + 1
#	comp<-data.exp[,c(index:(index+4))]
#	pairs(comp, upper.panel=panel.cor, cex.labels=1.0)
#}

# Hierarchical clustering
# set sample color
#groups<-tpmgroup
#group.sum<-summary(factor(groups, levels=unique(groups)))
#basecolor<-brewer.pal(length(unique(groups)), name="Accent")
#samplecolors=vector()
#for ( i in 1:length(basecolor))
#{
#    samplecolors=c(samplecolors, rep(basecolor[i], each=unname(group.sum[i])))
#}
#hc.title<-c("Hierarchical Clustering of All Samples")
#data.dist <- function(x) dist(x, method="euclidean")
#data.hclust <- function(d) hclust(d, method="ward")
#hmcol=colorRampPalette(brewer.pal(9, "RdBu" )) (100)
#nc<-length(names(data.exp))
#data.hm <-heatmap.2( as.matrix(log2(data.exp+1)),trace="none", distfun=data.dist,  hclustfun=data.hclust, cexCol = 0.2 + 1/log10(nc), 
#             margin=c(10, 10), main=hc.title, ColSideColors=samplecolors, labRow="", labCol=NULL, scale="row", col=rev(hmcol) )

# PCA for All SampleGroup
pca.title<-c("PCA of All Samplegroups")

cohort<-tpmgroup
#cohort<-tpmmodel
pca <- prcomp(t(log2(data.exp+1)))
score<-data.frame(cohort, pca$x[,1:3])
#colorbycohort<-factor(cohort)
colorbySamplegroups<-factor(cohort)
label<-tpmname
#ggplot(data=score, aes(x=PC1, y=PC2, colour=colorbycohort, size=4)) + 
ggplot(data=score, aes(x=PC1, y=PC2, colour=colorbySamplegroups, size=4)) + 
	geom_point(alpha = 1) +
#	xlim(-170, 100) +
#	geom_text(aes(label=label), fontface="bold", hjust=1,vjust=0, size=3, position=position_jitter(width=2,height=1) ) + 
	theme_bw() +
	theme(axis.text=element_text(face="bold"),
		  axis.title=element_text(face="bold"),
		  panel.border=element_rect(colour = "black", fill=NA, size=2))+
    scale_size_identity(guide="none") +
    ggtitle(pca.title)




# PCA for All Cellline
pca.title<-c("PCA of All Celllines")

#cohort<-tpmgroup
cohort<-tpmmodel
pca <- prcomp(t(log2(data.exp+1)))
score<-data.frame(cohort, pca$x[,1:3])
#colorbycohort<-factor(cohort)
colorbyCelllines<-factor(cohort)
label<-tpmname
#ggplot(data=score, aes(x=PC1, y=PC2, colour=colorbycohort, size=4)) + 
pcaplot<-ggplot(data=score, aes(x=PC1, y=PC2, colour=colorbyCelllines, size=4)) + 
	geom_point(alpha = 1) +
#	xlim(-170, 100) +
#	geom_text(aes(label=label), fontface="bold", hjust=1,vjust=0, size=3, position=position_jitter(width=2,height=1) ) + 
	theme_bw() +
	theme(axis.text=element_text(face="bold"),
		  axis.title=element_text(face="bold"),
		  panel.border=element_rect(colour = "black", fill=NA, size=2))+
    scale_size_identity(guide="none") +
    ggtitle(pca.title)

print(pcaplot)






if (length(unique(sampleinfo$CellLine)) > 1 )
{
	totalcell <- unique(sampleinfo$CellLine) 
	for (i in 1:length(totalcell) )
	{
		cell <- totalcell[i]	
		pca.title<-cell
		#fastqid<-sampleinfo$FastqID[sampleinfo$CellLine %in% cell]
                fastqid<-sampleinfo$SampleName[sampleinfo$CellLine %in% cell]
		graph.data<-data.exp[,c(fastqid)]
		tpmgroup<-samplegroup[match(fastqid, sampleid)]
		tpmname<-samplename[match(fastqid, sampleid)]
		cohort<-tpmgroup
		pca <- prcomp(t(log2(graph.data+1)))
		score<-data.frame(cohort, pca$x[,1:3])
		#colorbycohort<-factor(cohort)
		colorwithinCellline<-factor(cohort)
		label<-tpmname
		#pcaplot<-ggplot(data=score, aes(x=PC1, y=PC2, colour=colorbycohort, size=4)) + 
		pcaplot<-ggplot(data=score, aes(x=PC1, y=PC2, colour=colorwithinCellline, size=4)) + 
			   geom_point(alpha = 1) +
			#	xlim(-170, 100) +
			#	geom_text(aes(label=label), fontface="bold", hjust=1,vjust=0, size=3, position=position_jitter(width=2,height=1) ) + 
				theme_bw() +
				theme(axis.text=element_text(face="bold"),
		  		axis.title=element_text(face="bold"),
		  		panel.border=element_rect(colour = "black", fill=NA, size=2))+
    			scale_size_identity(guide="none") +
    			ggtitle(pca.title)
		print(pcaplot)
	}
}

dev.off()
