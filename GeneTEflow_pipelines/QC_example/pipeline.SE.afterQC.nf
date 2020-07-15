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
    
    trimmedfastq3=Channel
    .fromFilePairs( params.trimmedreads,size:-1 )
    .ifEmpty { exit 1, error "Cannot find any reads matching: ${params.trimmedreads}" }
} 
else if (params.adapter_trim_tag == "N" )  {
    trimmedfastq3=Channel
    .fromFilePairs( params.reads,size:-1 )
    .ifEmpty { exit 1, error "Cannot find any reads matching: ${params.reads}" }

}
else {
    println("Please set adapter trimming tag!")
    exit 1
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

genecounts1= Channel.fromPath(params.genecounts)


if (params.DESeq_run_tag == "Y" ) {
    if(params.DESeq_replicates == "N"){
       params.DESeq2_RSEM_cmd = "/RNASeq/scripts/RSEM_STAR.tximport.DESeq2.without_replicates.R"       
    }

    process DESeq2FromRSEM{
       publishDir params.AllResultsDir, mode: 'copy'
       publishDir params.ReportDir, mode: 'copy'


       input:
       file genecounts from genecounts1.collect()
       file sampleinfotxt from sampleinfo5
       file samplecomparetxt from samplecompare1

       output:
       file "*deseq2.RSEM_STAR.out.txt" into DESeq2_RSEM_result1,DESeq2_RSEM_result2

       """
       Rscript $params.DESeq2_RSEM_cmd  $sampleinfotxt  $samplecomparetxt  $params.compare_group  $params.tx2genefile  > DESeq2FromRSEM.log 2>DESeq2FromRSEM.err.log
       """
    }
}



GeneTPMRData2 = file(params.tpmdata)

if ((params.DESeq_run_tag == "Y" ) && (params.analysis_pipeline_run_tag == "Y"))
{
    process sigGene{
         publishDir params.ReportDir, mode: 'copy'
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file sampleinfotxt from sampleinfo6
         file samplecomparetxt from samplecompare2
         file deseqout from DESeq2_RSEM_result1
         file tpmdata from GeneTPMRData2

         output:
         file "*.siggenelist" into siggenelist
         file "*.txt" into outputtxt
         file "*.pdf" into outputpdf

         """    
         Rscript $params.sigGene_cmd $deseqout  $tpmdata  $sampleinfotxt  $samplecomparetxt  $params.deseq.fdr.TE  $params.deseq.log2FC.TE  $params.deseq.gmean.TE  >    sigGene.log   2>&1
         """

    }
}


if ((params.DESeq_run_tag == "Y" ) && (params.analysis_pipeline_run_tag == "Y")){

    process genGSEAinput{
         publishDir params.ReportDir, mode: 'copy'
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file samplecomparetxt from samplecompare3
         file deseqout from DESeq2_RSEM_result2

         output:
         file "samplecompare.gsea.txt" into samplecompareGSEAtxt
         file "*.rnk" into rnkfiles
 
         """
         Rscript $params.genGSEAinput_cmd $deseqout   $samplecomparetxt   >   genGSEAinput.log   2>&1
         """

     }
}



if ((params.DESeq_run_tag == "Y" ) && (params.analysis_pipeline_run_tag == "Y")){

    process submitGSEA{
         publishDir params.ReportDir, mode: 'copy'
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file rnkfile from rnkfiles.flatten()

         output:
         file "*.zip" into zipGSEAPrerankedfiles
  
         script:
         rptlabel = rnkfile[0].toString() - ~/.rnk/
    

         """
         java -cp $params.GSEA_jar -Xmx1024m xtools.gsea.GseaPreranked  -gmx $params.GSEA_gmx -norm meandiv -nperm 1000 -rnk $rnkfile  -scoring_scheme classic  -rpt_label $rptlabel  -create_svgs false -make_sets true -plot_top_x 40  -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report true -out ./ -gui false    >   genGSEAinput.log   2>&1

         zip -r ${rptlabel}.GseaPreranked.zip  ${rptlabel}.GseaPreranked.*
         """

     }
}





if (params.TE_pipeline_run_tag == "Y" ) {
    process squireFetch  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:

         output:
         file "squire_fetch" into squire_fetch

         """
         squire Fetch  -c -b  $params.squireFetch.genome  -f -r -g -x -p $params.cpu  -v >squire_Fetch.log   2>squire_Fetch.err.log
         """

     }

     process squireClean  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file squire_fetch

         output:
         file "squire_clean" into squire_clean

         """
         squire Clean  -b  $params.squireFetch.genome   -v >squire_Clean.log   2>squire_Clean.err.log
         """

     }


     process squireMap  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file squire_fetch
         file squire_clean
         set pair_id, file(reads) from trimmedfastq3

         output:
         set pair_id, "*.squire_map_folder" into squire_map

         """
         squire  Map  -1   ${pair_id}_1.clipped.fastq.gz      -o  ${pair_id}.squire_map_folder  -r  150  -n  ${pair_id}  -b  $params.squireFetch.genome   --gtf  squire_fetch/hg38_refGene.gtf  -p  $params.cpu  -v  >squire_Map_${pair_id}.log  2>squire_Map_${pair_id}.err.log
         """

     }


     process squireCountProcess  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file squire_fetch
         file squire_clean
         set pair_id, file(squire_map) from squire_map

         output:
         set pair_id, "*.squire_count_folder" into squire_count
         set pair_id, "*.TE_countTable.txt" into squire_TE1,squire_TE2

         """
         squire  Count  -m ${pair_id}.squire_map_folder   -o  ${pair_id}.squire_count_folder   -r  150 -n  ${pair_id}   -b  $params.squireFetch.genome   -s 1  -p  $params.cpu  -v>squire_Count_${pair_id}.log  2>squire_Count_${pair_id}.err.log
         
         cat ${pair_id}.squire_count_folder/${pair_id}_TEcounts.txt|cut -f 4,16 >  ${pair_id}.TE_countTable.txt
         """

     }

     process squireTE_list  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file TE_countTable from squire_TE1.collect()

         output:
         file "TE_id.list" into squire_TE_list

         """
         cat *.TE_countTable.txt |cut -f 1|grep -v "TE_ID"|sort|uniq >  TE_id.list
         """

     }

     process squire_TE_CountTable  {
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file squire_TE_list from squire_TE_list
         set pair_id, file(squire_count) from squire_TE2

         output:
         set pair_id, "*.TEcount.txt" into squire_TEcountTable

         """
         perl  /RNASeq/scripts/parse_list_from_TEcounts.pl  ${squire_TE_list}  ${pair_id}.TE_countTable.txt  > ${pair_id}.TEcount.txt
         """

     }



process genTEMatrixCounts {
    publishDir params.AllResultsDir, mode: 'copy'
    publishDir params.ReportDir, mode: 'copy'

    input:
    file TEcounts from squire_TEcountTable.collect()
    file sampleinfotxt from sampleinfo10

    output:
    file 'all.sample.Counts.TE.results.txt' into TECounts

    """
    perl /RNASeq/scripts/genTECountsMatrix.pl    .  $sampleinfotxt  counts > genTEMatrixCounts.log 2>genTEMatrixCounts.err.log
    """
}


process genTERData {
    publishDir params.AllResultsDir, mode: 'copy'

    input:
    file sampleinfotxt from sampleinfo11
    file TEcounts from TECounts

    output:
    file 'all.sample.TE.counts.rda' into TECountsRData
    file 'all.sample.TE.counts.txt' into TECountstxt
    file 'all.sample.TE.tpm.rda' into TETPMRData
    file 'all.sample.TE.tpm.txt' into TETPMtxt



    """
    Rscript /RNASeq/scripts/genTERDA.R  $sampleinfotxt  $TEcounts  >genRData.log  2>genRData.err.log
    mv all.sample.counts.rda all.sample.TE.counts.rda
    mv all.sample.counts.txt  all.sample.TE.counts.txt
    mv all.sample.tpm.rda all.sample.TE.tpm.rda
    mv all.sample.tpm.txt  all.sample.TE.tpm.txt
    """
}


    process DESeq2TE{
       publishDir params.AllResultsDir, mode: 'copy'
       publishDir params.ReportDir, mode: 'copy'


       input:
       file TEcounts from TECountsRData
       file sampleinfotxt from sampleinfo12
       file samplecomparetxt from samplecompare5

       output:
       file "*deseq2.TE.out.txt" into DESeq2_TE_result1

       """
       Rscript /RNASeq/scripts/DESeq2.TE.findDE  $TEcounts  $sampleinfotxt  $samplecomparetxt  $params.compare_group   > DESeq2TE.log 2>DESeq2TE.err.log
       """
    }

    process sigTE{
         publishDir params.ReportDir, mode: 'copy'
         publishDir params.AllResultsDir, mode: 'copy'

         input:
         file sampleinfotxt from sampleinfo13
         file samplecomparetxt from samplecompare6
         file deseqout from DESeq2_TE_result1
         file tpmdata from TETPMRData

         output:
         file "*.sigTElist" into siggenelistTE
         file "*.txt" into outputtxtTE
         file "*.pdf" into outputpdfTE

         """
         Rscript /RNASeq/scripts/sigTE.R   $deseqout  $tpmdata  $sampleinfotxt  $samplecomparetxt  0.05  0.6  5  >    sigTE.log   2>&1
         """

    }


}







workflow.onComplete { 
        
        endtime = new java.util.Date()
        
	println ( workflow.success ? "\n\nRNA-Seq pipeline is done at $endtime\n\n" : "Oops .. something went wrong!!!" )
}


