#!/bin/Rscript
#
# This script generates summary information of statistically regulated genes
# from DESeq2 analysis
# Author: Wenyan Zhong
# last revision date: 07/09/2018
# run it in the */results directory  
#
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(tables)
library(gtable)
library("RColorBrewer")
library("gplots")
library("ComplexHeatmap")
library(circlize)
library(VennDiagram)

## input files
args = commandArgs(TRUE)

#sampleinfo.file <- args[1]   #sample information txt file sampleinfo.txt
#tpm.file <- args[2]         #TPM file all.sample.TPM.genes.results

#allDE.file<-c("all.deseq2.out.txt")  # deseq2 output file
#tpm.file<-c("all.sample.tpm.rda")	# Rdata of all tpm data
#sampleinfo.file<-c("../sampleinfo/sampleinfo.txt")  # sample information file
#samplecompare.file<-c("../sampleinfo/samplecompare.txt") # sample comparison file
#fdr.cutoff<-0.05 	# cut off for fdr
#ratio.cutoff<-1	    # cut off for log2FC 
#mean.cutoff<-50		# cut off for max counts of group mean 


allDE.file <- args[1]                       # deseq2 output file
tpm.file <- args[2]                         # Rdata of all tpm data
sampleinfo.file <- args[3]                  # sample information file
samplecompare.file <- args[4]               # sample comparison file
fdr.cutoff <- as.double(args[5])           # cut off for fdr
print("FDR cutoff:")
print(fdr.cutoff)
ratio.cutoff <- as.double(args[6])         # cut off for log2FC
print("log2FC cutoff:")
print(ratio.cutoff)
mean.cutoff <- as.double(args[7])          # cut off for max counts of group mean
print("max counts of group mean cutoff:")
print(mean.cutoff)


# read in differential data
allDE<-read.table(file=allDE.file, header=TRUE, sep="\t", check.names=F, stringsAsFactor=F)
# get sample info
sampleinfo<-read.table(sampleinfo.file, sep="\t", check.names=F, stringsAsFactor=F, header=TRUE )
samplename<-sampleinfo$SampleName
fastqid<-sampleinfo$FastqID
# get sample compare info
samplecompare<-read.table(samplecompare.file, header=T, sep="\t", check.names=F, stringsAsFactor=F)
# read in tpm data
data.exp<-readRDS(tpm.file)
tpmid<-names(data.exp)
#tpmname<-samplename[match(tpmid, fastqid)]
#names(data.exp)<-tpmname

## extract data for generating significant gene list
samplecompare<-samplecompare[samplecompare$Stats %in% "y",]
compareGroups<-unique(c(samplecompare$Group1, samplecompare$Group2))
sampleinfo<-sampleinfo[sampleinfo$SampleGroup %in% compareGroups,]
compareSampleNames<-sampleinfo$SampleName
data.exp<-data.exp[, compareSampleNames]

# get sig genes for each comparisons
totalcompare<-nrow(samplecompare)
compare<-paste(samplecompare$Group2, samplecompare$Group1, sep=".vs.")
sig.sum<-list()
#siggenes<-vector()
siggeneslist<-list()
for ( i in 1:totalcompare )
{
	group1.mean.col<-paste(samplecompare$Group1[i], "Mean", sep=".")
	group2.mean.col<-paste(samplecompare$Group2[i], "Mean", sep=".")
	if ( any(j <- grep("^\\d+", group1.mean.col ) ))
	{
		group1.mean.col<-paste("X", group1.mean.col, sep="" )
	}
	
	if ( any(j <- grep("^\\d+", group2.mean.col ) ))
	{
		group2.mean.col<-paste("X", group2.mean.col, sep="")
	}
	group1.mean<-allDE[,group1.mean.col] 
	group2.mean<-allDE[,group2.mean.col]
	group.mean<-data.frame(group1.mean=group1.mean, group2.mean=group2.mean)
	group.mean$maxmean<-apply(group.mean, 1, max)
	ratio.i<-paste(compare[i], "Log2FC", sep=".")
	pval.i<-paste(compare[i], "pval", sep=".")
	fdr.i<-paste(compare[i], "padj", sep=".")
	sig.f<-!is.na(allDE[,ratio.i]) & !is.na(allDE[,fdr.i]) & allDE[,fdr.i]<=fdr.cutoff & (allDE[,ratio.i]>=ratio.cutoff | allDE[,ratio.i]<=-ratio.cutoff) & group.mean$maxmean >=mean.cutoff
	sig<-allDE[sig.f,]
	sig.up<-sig[sig[,ratio.i] >= ratio.cutoff,]
	sig.down<-sig[sig[,ratio.i] <= -ratio.cutoff,]
	sig.i.sum<-data.frame("UP"=nrow(sig.up), "DOWN"=nrow(sig.down), "Total"=nrow(sig) )
	sig.sum<-rbind(sig.sum, sig.i.sum)
	siggene.i<-sig$GeneSymbol
	siggeneslist[[compare[i]]]<-siggene.i
#	siggenes<-union(siggenes, siggene.i)
}
row.names(sig.sum)<-compare

## output summary of significant genes
output.sigsum.file<-c("all.sample.siggene.sum.txt")
write.table(file=output.sigsum.file, cbind(data.frame("Comparison"=row.names(sig.sum)), sig.sum), row.names=F, quote=F, sep="\t")

#output sig gene list
for ( k in 1:length(siggeneslist))
{
	genelist.filename<-paste(names(siggeneslist)[k], ".siggenelist", sep="")
	genelist<-data.frame(compname=siggeneslist[[k]])
	names(genelist)<-names(siggeneslist)[k]
	write.table(file=genelist.filename, genelist, sep="\t", row.names=F, quote=F)
}

# generates significant gene comparison graph for subset of samplegroups
compareSubGroup<-unique(samplecompare$SubGroup)

for ( i in 1:length(compareSubGroup))
{
	subsetGroup<-compareSubGroup[i]

	comparesubset<-samplecompare[samplecompare$SubGroup %in%  subsetGroup,]
	subset<-paste(comparesubset$Group2, ".vs.", comparesubset$Group1, sep="") 
	# get subset of sig genes
	siggenesubset=vector()
	for (j in 1:length(subset))
	{
		sig.j<-siggeneslist[[subset[j]]]
		siggenesubset<-union(siggenesubset, sig.j)
	}
	# get data for union of subset compare genes and samples in the subsets
	comparegroup<-unique(c(comparesubset$Group2, comparesubset$Group1))
	comparesampleinfo<-sampleinfo[sampleinfo$SampleGroup %in% comparegroup,]
	comparesample<-comparesampleinfo$SampleName
	data.exp.subset<-data.exp[row.names(data.exp) %in% siggenesubset, c(comparesample)]

	outputunion.sigfile<-paste("all.sig.gene.subset.", subsetGroup, ".tpm.txt", sep="")
	write.table(cbind("GeneSymbol"=row.names(data.exp.subset), data.exp.subset), outputunion.sigfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

	## generate graphs 
	grpahfile<-paste("all.sig.gene.subset",subsetGroup, fdr.cutoff, ratio.cutoff, "sum.pdf", sep="." )
	pdf(file=grpahfile)
	table.theme <- ttheme_default( 
					 core = list(fg_params=list(fontface="bold", fontsize=12,
								 parse=FALSE),
								 bg_params=list(fill=c("white", "lightblue"))),
					 colhead = list(fg_params=list( fontface="bold.italic",fontsize=14,
								 parse=FALSE),
								 bg_params=list(fill="lightpink"))
	)
	sig.title<-c("Number of Significantly \n Regulated Genes \n")
        #sig.title<-c("Number of DEG")	
	table <- tableGrob(sig.sum[subset,],theme=table.theme )
	#title <- textGrob(sig.title, y=unit(0.9,"npc")  , 
	title <- textGrob(sig.title, gp=gpar(fontsize=14,fontface="bold" ))

	sigcond<-paste("FDR<=", fdr.cutoff, "\n|logRatio|>=", ratio.cutoff, sep="")
	footnote <- textGrob(sigcond, x=0, hjust=0,
						 gp=gpar( fontface="italic"))

	padding <- unit(0.5,"line")
	table <- gtable_add_rows(table,  
							 heights = grobHeight(title) + padding,
							 pos = 0)
	table <- gtable_add_rows(table, 
							 heights = grobHeight(footnote)+ padding)
	table <- gtable_add_grob(table, list(title, footnote),
							 t=c(1, nrow(table)), l=c(1,2), 
							 r=ncol(table))

	grid.newpage()
	grid.draw(table)

	## generate venndiagram 
	numcolor=length(subset)
	if (numcolor <= 5 )
	{     
           if(numcolor < 3){fillcolor<-brewer.pal(3, "Set1")}
           else{fillcolor<-brewer.pal(numcolor, "Set1")}
        
           fillcolor<-fillcolor[1:numcolor]
		#fillcolor<-brewer.pal(numcolor, "Set1")
		cat.dist = rep(0.04, numcolor)
		venn<-venn.diagram(
		x = siggeneslist[subset],
		filename = NULL,
		height=600,
		width=600,
		lwd =1,
		lty=1,
		col = "black",
		fill = fillcolor,
		alpha = 0.5,
		label.col = c("white"),
		cex = 1.5,
		fontface = "bold",
		cat.col = fillcolor,
		cat.cex = 1.0,
		cat.fontfamily = "serif",
		cat.dist = cat.dist,
	#    cat.pos = c(-20, 20 ),
		euler.d = TRUE,
		scaled = TRUE,
		margin=0.05,
		main.cex=1.5,
		)

		grid.newpage()
		grid.draw(venn)
	}
	## generate HC Plots for significantly changed genes

	mat<-log2(data.exp.subset +1 )
	rs=rowSums(data.exp.subset)
	use =  (rs > 10)
	mat.rowscaled<-t(apply(data.exp.subset, 1, scale) )
	colnames(mat.rowscaled)<-colnames(data.exp.subset)
#	colnames(mat.rowscaled)<-gsub("_L(.*)$", "", colnames(mat.rowscaled), perl=T)
	ht<-Heatmap(mat.rowscaled, col = colorRamp2(c(-2,0,2), c("blue", "white", "red")),
            show_column_names =T,
			show_row_names = F,
			column_names_max_height = unit(7, "cm"),
            column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    #        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
    #        row_names_max_width = unit(16, "cm"),
            name="expression",
            heatmap_legend_param = list(
                                color_bar = "continuous",
                             #   legend_direction = "horizontal",
				legend_direction = "vertical",
                                legend_width = unit(10, "cm"),
                                #title_position = "lefttop")
				title_position = "topcenter")
    )
	plot.heatmap<-draw(ht, heatmap_legend_side = "right")
	dev.off()
}
