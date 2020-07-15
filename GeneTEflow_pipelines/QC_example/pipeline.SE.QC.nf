//Global declarations
////////////////////////////////
//Scripts:
params.fastqc = "fastqc"
params.rsem = "rsem-calculate-expression"
//params.RSEM_STAR_Index = "/RNASeq/ngsdb/RSEM_STAR_Index_hg38_UCSC/GRChg38_UCSC"
//params.DESeq2_RSEM_cmd = "/RNASeq/scripts/RSEM_STAR.tximport.DESeq2.without_replicates.R"
params.DESeq2_RSEM_cmd = "/RNASeq/scripts/RSEM_STAR.tximport.DESeq2.R"
params.tx2genefile = "/RNASeq/ngsdb/tximport_library/tx2gene.GRChg38_UCSC.csv"
params.trimmomatic = "trimmomatic"
params.adapterfa = "TruSeq3-SE"
params.genSampleInfo_cmd = "/RNASeq/scripts/genSampleInfo.R"
params.genGeneEXPMatrix_cmd="/RNASeq/scripts/genGeneEXPMatrix.pl"
//params.mapstat_rsem_cmd="/RNASeq/scripts/getmapstat.rsem.pl"
params.mapstat_star_cmd="/RNASeq/scripts/getmapstat.star.pl"
params.genRDA_cmd="/RNASeq/scripts/genRDA.R"
params.QC_cmd="/RNASeq/scripts/QC.R"
params.sigGene_cmd="/RNASeq/scripts/sigGene.R"
params.genGSEAinput_cmd="/RNASeq/scripts/genGSEAinput.R"
params.GSEA_jar="/RNASeq/apps/GSEA-3.0/gsea-3.0.jar"
params.GSEA_gmx="/RNASeq/ngsdb/GSEA_symbols/c2.all.v6.1.symbols.gmt"
params.GSEA_chipfile="/RNASeq/ngsdb/GSEA_symbols/GENE_SYMBOL.msigdb_v6.1.chip"
params.genGSEAExpinput_cmd="/RNASeq/scripts/genGSEAinput.exp.R"
///////////////////////////////




starttime = new java.util.Date()
println "\n\n\nRNA-Seq pipeline is starting at "+starttime+"\n\n\n";


if (params.adapter_trim_tag == "Y" ) {
    
    read_pairs=Channel
    .fromFilePairs( params.reads ,size:-1 )
    .ifEmpty { exit 1, error "Cannot find any reads matching: ${params.reads}" }


    process Trim_Adapters {
        publishDir params.TrimmedDir, mode: 'copy'

        input:
        set pair_id, file(reads) from read_pairs

        output:
        set pair_id, "*.clipped.fastq.gz"  into trimmedfastq1, trimmedfastq2, trimmedfastq3

        script:
        """
        $params.trimmomatic  SE -threads $params.cpu -phred33  ${pair_id}_R1${params.reads_suffix}    ${pair_id}_1.clipped.fastq.gz   ILLUMINACLIP:$params.adapterfa:2:30:10 > trimmomatic.log 2>&1
        """
    }
} 
else if (params.adapter_trim_tag == "N" )  {
    trimmedfastq1=Channel
    .fromFilePairs( params.reads,size:-1 )
    .ifEmpty { exit 1, error "Cannot find any reads matching: ${params.reads}" }

    trimmedfastq2=Channel
    .fromFilePairs( params.reads ,size:-1)
    .ifEmpty { exit 1, error "Cannot find any reads matching: ${params.reads}" }

}
else {
    println("Please set adapter trimming tag!")
    exit 1
}




process FASTQC {
    publishDir params.FastQCdir, mode: 'copy'

    input:
    set pair_id, file(reads) from trimmedfastq1

    output:
    set pair_id, "*.zip" into fastqczip
    set pair_id, "*.html" into fastqchtml


    """    
    $params.fastqc  --noextract -f  fastq    -o   ./   ${pair_id}_1.clipped.fastq.gz     >    ${pair_id}_1.fastqc.log   2>&1
    """

}


fasta_ref = file(params.fasta)
gencode_annotation_gtf = file(params.gtf)


process Make_STARindex {
    publishDir params.AllResultsDir, mode: 'copy'

    input:
    file fasta_ref
    file gencode_annotation_gtf
            
    output:
    file "index_folder" into star_index
           
    """
    mkdir index_folder
    rsem-prepare-reference  -p $params.cpu   --gtf  ${gencode_annotation_gtf}     --star  ${fasta_ref}   index_folder/GRChg38_UCSC  > RSEM_STAR_Index_hg38_UCSC.log  2>RSEM_STAR_Index_hg38_UCSC.err.log
    """
}



process RSEMwithSTAR {
    publishDir params.AllResultsDir, mode: 'copy'
    
    input:
    set pair_id, file(reads) from trimmedfastq2
    file RSEM_STAR_Index from star_index
 
    output:
    set pair_id, "*.genes.results" into RSEMresult1,RSEMresult2,RSEMresult3,RSEMresult4
    set pair_id, "*.isoforms.results" into RSEMresult_isoform
    set pair_id, "*.cnt" into RSEMresult_cnt
    set pair_id, "*.RSEM_Output.STAR.log" into RSEMresult_STAR_log
 
    """    
    $params.rsem  -p  $params.cpu  --forward-prob $params.forwardprob  --sort-bam-by-coordinate  --keep-intermediate-files    --star    --star-gzipped-read-file    ${pair_id}_1.clipped.fastq.gz      index_folder/GRChg38_UCSC  ${pair_id}.RSEM_Output  >  ${pair_id}.RSEM.log  2> ${pair_id}.RSEM.err.log
 
    cp ${pair_id}.RSEM_Output.stat/${pair_id}.RSEM_Output.cnt  ${pair_id}.RSEM_Output.cnt
    cp ${pair_id}.RSEM_Output.temp/${pair_id}.RSEM_OutputLog.final.out  ${pair_id}.RSEM_Output.STAR.log
    """

}


sampleinfoxlsx = file(params.sampleinfoxlsx)

process genSampleInfo {
    publishDir params.sampleinfoDir, mode: 'copy'
    
    input:
    file sampleinfoxlsx

    output:
    file 'sampleinfo.txt' into sampleinfo1,sampleinfo2,sampleinfo3,sampleinfo4,sampleinfo5,sampleinfo6,sampleinfo7,sampleinfo8,sampleinfo9,sampleinfo10,sampleinfo11,sampleinfo12,sampleinfo13
    file 'samplecompare.txt' into samplecompare1,samplecompare2,samplecompare3,samplecompare4,samplecompare5,samplecompare6

    """
    Rscript $params.genSampleInfo_cmd  $sampleinfoxlsx  $params.file.extention  $params.sample.manifest.sheetname  $params.samplecompare.sheetname  > genSampleInfo.log 2>genSampleInfo.err.log
    """
}


process genGeneEXPMatrixTPM {
    publishDir params.AllResultsDir, mode: 'copy'
    publishDir params.ReportDir, mode: 'copy'

    input:
    file genecounts from RSEMresult2.collect()
    file sampleinfotxt from sampleinfo1

    output:
    file 'all.sample.TPM.genes.results.txt' into GeneTPM

    """
    perl $params.genGeneEXPMatrix_cmd  .  $sampleinfotxt  tpm > genGeneEXPMatrixTPM.log 2>genGeneEXPMatrixTPM.err.log
    """
}




process genGeneEXPMatrixCounts {
    publishDir params.AllResultsDir, mode: 'copy'
    publishDir params.ReportDir, mode: 'copy'

    input:
    file genecounts from RSEMresult4.collect()
    file sampleinfotxt from sampleinfo7

    output:
    file 'all.sample.Counts.genes.results.txt' into GeneCounts

    """
    perl $params.genGeneEXPMatrix_cmd  .  $sampleinfotxt  counts > genGeneEXPMatrixCounts.log 2>genGeneEXPMatrixCounts.err.log
    """
}







process getmapstatFromSTAR {
    publishDir params.AllResultsDir, mode: 'copy'
    publishDir params.ReportDir, mode: 'copy'

    input:
    file mapcountslog from RSEMresult_STAR_log.collect()
    file sampleinfotxt from sampleinfo8

    output:
    file 'all.sample.star_map.stat.txt' into mapstatFromSTAR

    """
    perl $params.mapstat_star_cmd  .  $sampleinfotxt   > all.sample.star_map.stat.txt  2>getmapstat.err.log
    """
}


process genRData {
    publishDir params.AllResultsDir, mode: 'copy'

    input:
    file genetpm from GeneTPM
    file sampleinfotxt from sampleinfo3
    file genecounts from GeneCounts

    output:
    file 'all.sample.tpm.rda' into GeneTPMRData1,GeneTPMRData2,GeneTPMRData3
    file 'all.sample.tpm.txt' into GeneTPMtxt
    file 'all.sample.counts.rda' into GeneCountsRData
    file 'all.sample.counts.txt' into GeneCountstxt
    

    """
    Rscript $params.genRDA_cmd  $sampleinfotxt  $genetpm  $genecounts  >genRData.log  2>genRData.err.log
    """
}



process QC {
    publishDir params.AllResultsDir, mode: 'copy'
    publishDir params.ReportDir, mode: 'copy'


    input:
    file genetpmrda from GeneTPMRData1
    file sampleinfotxt from sampleinfo4

    output:
    file 'all.sample.TPM.stat.txt' into GeneTPMSTAT
    file 'Sample.QC.pdf' into GeneQCpdf

    """
    Rscript $params.QC_cmd  $sampleinfotxt  $genetpmrda  >QC.log  2>QC.err.log
    """
}








workflow.onComplete { 
        
        endtime = new java.util.Date()
        
	println ( workflow.success ? "\n\nRNA-Seq pipeline is done at $endtime\n\n" : "Oops .. something went wrong!!!" )
}


