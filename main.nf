#!/usr/bin/env nextflow 

/*
========================================================================================
                            m6APipe
========================================================================================
 * m6Apipe
 * Homepage / Documentation

 */

/*
 * to be added
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 * Zhu Kaiyu <zhuky5@mail2.sysu.edu.cn>:
 */

// requirement:
// - fastp/fastqc
// - STAR/tophat2/bowtie2/hisat2
// - samtools/rseqc
// - MeTPeak/Macs2/MATK/meyer
// - RobustRankAggreg/bedtools
// - Cufflinks/DESeq2/EdgeR
// - MeTDiff/QNB/MATK
// - Homer/DREME
// - SRAMP/MATK
// - igvtools

//pre-defined functions for render command
//=======================================================================================
def helpMessage() {
    println LikeletUtils.sysucc_ascii()
    log.info"""
    =========================================
     nf-core/m6Apipe v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow path/to/m6APipe/main.nf --readPaths './data/' --designfile -profile standard,docker 

    Mandatory arguments:
      --readPaths                   Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      --designfile                  format:filename,control_or_treated,ip_or_input,tag_id
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test
    
    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF reference

    Options:
      --inputformat                 fastq.gz;fastq default = fastq
      --singleEnd                   Specifies that the input is single end reads
      --tophat2_index               Path to tophat2 index, eg. "path/to/Tophat2Index/*"
      --hisat2_index                Path to hisat2 index, eg. "path/to/Hisat2Index/*"
      --bwa_index                   Path to bwa index, eg. "path/to/BwaIndex/*"
      --star_index                  Path to star index, eg. "path/to/StarIndex/"
      --skip_qc                     Skip all QC steps                        
      --skip_expression             Skip all differential expression analysis steps
      --skip_peakCalling            Skip all Peak Calling steps
      --skip_diffpeakCalling        Skip all Differential methylation analysis

    
    Other options:
      --outdir                      The output directory where the results will be saved, defalut = $baseDir/results
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_genebody_coverage      Skip calculating genebody coverage  
      --skip_cufflinks              Skip the cufflinks process of differential expression analysis steps
      --skip_edger                  Skip the EdgeR process of differential expression analysis steps
      --skip_deseq2                 Skip the deseq2 process of differential expression analysis steps
      --skip_metpeak                Skip the metpeak process of Peak Calling steps
      --skip_macs2                  Skip the macs2 process of Peak Calling steps
      --skip_matk                   Skip the matk process of Peak Calling steps
      --skip_metdiff                Skip the metdiff process of Differential methylation analysis
      --skip_QNB                    Skip the QNB process of Differential methylation analysis
      --skip_diffmatk               Skip the matk process of Differential methylation analysis
 

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.call = false
params.email = false
params.plaintext_email = false
params.seqCenter = false
params.help = false


// Preset trimming options
/*params.pico = false
if (params.pico){
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}
*/
// Validate inputs

if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, LikeletUtils.print_red("Fasta file not found: ${params.fasta}")
}else {
    exit 1, LikeletUtils.print_red("No reference genome specified!")
}
if( params.gtf ){
    gtf = file ( params.gtf )
    if( !gtf.exists() ) exit 1, LikeletUtils.print_red("gtf not found: ${params.gtf}")
} else {
    exit 1, LikeletUtils.print_red("No GTF annotation specified!")
}
if( params.designfile ) {
    designfile = file(params.designfile)
    if( !designfile.exists() ) exit 1, LikeletUtils.print_red("Design file not found: ${params.designfile}")
}else{
    exit 1, LikeletUtils.print_red("No Design file specified!")
}
if( params.comparefile == "two_group" ){
    comparefile = false
    compareLines = Channel.from("two_group")
}else if( params.comparefile){
    comparefile = file(params.comparefile)
    if( !comparefile.exists() ) exit 1, print_red("Compare file not found: ${params.comparefile}")
    compareLines = Channel.from(comparefile.readLines())
}else {
    comparefile = false
    compareLines = Channel.from("")
}
compareLines.into{
    compareLines_for_DESeq2; compareLines_for_edgeR; compareLines_for_plot;
    compareLines_for_diffm6A; compareLines_for_arranged_result
}
compareLines_for_plot.subscribe{print it}
// Validate the params of skipping Aligners Tools Setting
if( params.aligners == "none" ){
    skip_aligners = true
    skip_bwa = true
    skip_tophat2 = true
    skip_hisat2 = true
    skip_star = true
}else if( params.aligners == "star" ){
    skip_aligners = false
    skip_bwa = true
    skip_tophat2 = true
    skip_hisat2 = true
    skip_star = false
}else if( params.aligners == "hisat2" ){
    skip_aligners = false
    skip_bwa = true
    skip_tophat2 = true
    skip_hisat2 = false
    skip_star = true
}else if( params.aligners == "tophat2" ){
    skip_aligners = false
    skip_bwa = true
    skip_tophat2 = false
    skip_hisat2 = true
    skip_star = true
}else if( params.aligners == "bwa" ){
    skip_aligners = false
    skip_bwa = false
    skip_tophat2 = true
    skip_hisat2 = true
    skip_star = true
}else{
    exit 1, LikeletUtils.print_red("Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'tophat2', 'bwa'")
}
if( params.expression_analysis_mode == "edgeR" ){
    params.skip_edger = false
    params.skip_deseq2 = true
    params.skip_cufflinks = true
    params.skip_expression = false
}else if( params.expression_analysis_mode == "DESeq2" ){
    params.skip_edger = true
    params.skip_deseq2 = false
    params.skip_cufflinks = true
    params.skip_expression = false
}else if( params.expression_analysis_mode == "Cufflinks" ){
    params.skip_edger = true
    params.skip_deseq2 = true
    params.skip_cufflinks = false
    params.skip_expression = false
}else if( params.expression_analysis_mode == "none" ){
    params.skip_edger = true
    params.skip_deseq2 = true
    params.skip_cufflinks = true
    params.skip_expression = true
}else{
    exit 1, LikeletUtils.print_red("Invalid expression_analysis_mode option: ${params.expression_analysis_mode}. Valid options: 'edgeR', 'DESeq2', 'none'")
}
/*
 * Create a channel for input read files
 */
if( params.readPaths && !skip_aligners ){
    if( params.singleEnd ){
        Channel
            .fromFilePairs( "${params.readPaths}/*.{fastq,fastq.gz}", size: 1 ) 
            .ifEmpty { exit 1, LikeletUtils.print_red("readPaths was empty - no fastq files supplied: ${params.readPaths}")}
            .into{ raw_data; raw_bam }
    }
    else if ( !params.singleEnd ){
        Channel
            .fromFilePairs( "${params.readPaths}/*{1,2}.{fastq,fastq.gz}", size: 2 ) 
            .ifEmpty { exit 1, LikeletUtils.print_red("readPaths was empty - no fastq files supplied: ${params.readPaths}") }
            .into{ raw_data; raw_bam }
    }
    else {
        exit 1, println LikeletUtils.print_red("The param 'singleEnd' was not defined!")
    }
}else if( params.readPaths && skip_aligners ){
    Channel
        .fromPath( "${params.readPaths}/*.bam") 
        .ifEmpty { exit 1, LikeletUtils.print_red("readPaths was empty - no bam files supplied: ${params.readPaths}")}
        .into{ raw_data; raw_bam }
} 
else{
    println LikeletUtils.print_red("readPaths was empty: ${params.readPaths}")
}
/*
========================================================================================
                         showing the process and files
========================================================================================
*/
println LikeletUtils.sysucc_ascii()
println LikeletUtils.print_purple("============You are running m6APipe with the following parameters===============")
println LikeletUtils.print_purple("Checking parameters ...")

println LikeletUtils.print_yellow("=====================================Reads types================================")
println (LikeletUtils.print_yellow("SingleEnd :                     ") + LikeletUtils.print_green(params.singleEnd))
println (LikeletUtils.print_yellow("Stranded :                      ") + LikeletUtils.print_green(params.stranded))
println (LikeletUtils.print_yellow("gzip :                          ") + LikeletUtils.print_green(params.gzip))

println LikeletUtils.print_yellow("====================================Mode selected==============================")
println (LikeletUtils.print_yellow("aligners :                      ") + LikeletUtils.print_green(params.aligners))
println (LikeletUtils.print_yellow("peakCalling_mode :              ") + LikeletUtils.print_green(params.peakCalling_mode))
println (LikeletUtils.print_yellow("peakMerged_mode :               ") + LikeletUtils.print_green(params.peakMerged_mode))
println (LikeletUtils.print_yellow("expression_analysis_mode :      ") + LikeletUtils.print_green(params.expression_analysis_mode))
println (LikeletUtils.print_yellow("methylation_analysis_mode :     ") + LikeletUtils.print_green(params.methylation_analysis_mode))

println LikeletUtils.print_yellow("==================================Input files selected==========================")
println (LikeletUtils.print_yellow("Reads Path :                    ") + LikeletUtils.print_green(params.readPaths))
println (LikeletUtils.print_yellow("fasta file :                    ") + LikeletUtils.print_green(params.fasta))
println (LikeletUtils.print_yellow("Gtf file :                      ") + LikeletUtils.print_green(params.gtf))
println (LikeletUtils.print_yellow("Design file :                   ") + LikeletUtils.print_green(params.designfile))
println (LikeletUtils.print_yellow("Compare file :                  ") + LikeletUtils.print_green(params.comparefile))

println LikeletUtils.print_yellow("==================================Skip model selected==========================")
println (LikeletUtils.print_yellow("Skip samtools sort :            ") + LikeletUtils.print_green(params.skip_sort))
println (LikeletUtils.print_yellow("Skip expression analysis :      ") + LikeletUtils.print_green(params.skip_expression))
println (LikeletUtils.print_yellow("Skip peakCalling :              ") + LikeletUtils.print_green(params.skip_peakCalling))
println (LikeletUtils.print_yellow("Skip diffpeakCalling :          ") + LikeletUtils.print_green(params.skip_diffpeakCalling))
println (LikeletUtils.print_yellow("Skip annotation :               ") + LikeletUtils.print_green(params.skip_annotation))
println (LikeletUtils.print_yellow("Skip qc :                       ") + LikeletUtils.print_green(params.skip_qc))

println LikeletUtils.print_yellow("==================================Output files directory========================")
println (LikeletUtils.print_yellow("Output directory :              ") + LikeletUtils.print_green(params.outdir))


/*
========================================================================================
                             check or build the index
========================================================================================
*/ 
/*
 * PREPROCESSING - Build BED12 file
 * NEED gtf.file
 */
process makeBED12 {
    label 'build_index'
    tag "gtf2bed12"
    publishDir path: { params.saveReference ? "${params.outdir}/Genome/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file gtf
    
    output:
    file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

    script:      
    """
    bash ${baseDir}/bin/gtf2bed12.sh $gtf
    """
}
/*
 * PREPROCESSING - Create chromsizesfile for meyer
 * NEED genome.fa
 */
if( params.chromsizesfile ) {
    chromsizesfile = file(params.chromsizesfile)
    if( !chromsizesfile.exists() ) exit 1, LikeletUtils.print_red("Chromsizes file not found: ${params.chromsizesfile}")
}else if( params.fasta ){
    process MeyerPrepration{
        label 'build_index'
        tag "MeyerPrepration"

        input:
        file fasta

        output:
        file "chromsizes.file" into chromsizesfile

        when:
        !params.chromsizesfile && !params.skip_meyer && !params.skip_peakCalling

        shell:
        '''
        cat !{fasta} | awk 'BEGIN{len=""}{if($0~">"){split($0,ID,"[> ]");printf len"ABC"ID[2]"\\t";len=0}else{len=len+length($0)}}END{print len}' |sed 's/ABC/\\n/g' |awk NF > chromsizes.file
        '''
    }
}else {
    exit 1, println LikeletUtils.print_red("Chromsizes file cannot build due to lack of ")
}
/*
 * PREPROCESSING - Build TOPHAT2 index
 * NEED genome.fa
 */
if( params.tophat2_index && !skip_tophat2 ){
    tophat2_index = Channel
        .fromPath( params.tophat2_index )
        .ifEmpty { exit 1, "Tophat2 index not found: ${params.tophat2_index}" }
}else if( params.fasta ){
    process MakeTophat2Index {
        label 'build_index'
        tag "tophat2_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/": params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        input:
        file fasta

        output:
        file "Tophat2Index/*" into tophat2_index

        when:
        !skip_tophat2 && !skip_aligners

        script:
        tophat2_index = "Tophat2Index/" + fasta.baseName.toString()
        """
        mkdir Tophat2Index
        bowtie2-build -f $params.fasta $tophat2_index
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no Tophat2 Index")
}

/*
 * PREPROCESSING - Build HISAT2 index
 * NEED genome.fa genes.gtf snp.txt/vcf
 */
if( params.hisat2_index && !skip_hisat2 ){
    hisat2_index = Channel
        .fromPath(params.hisat2_index)
        .ifEmpty { exit 1, "hisat2 index not found: ${params.hisat2_index}" }
}else if( params.fasta){
    process MakeHisat2Index {
        label 'build_index'
        tag "hisat2_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/ " : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'        
        input:
        file fasta
        file gtf
        //file params.snp

        output:
        file "Hisat2Index/*" into hisat2_index

        when:
        !skip_hisat2 && !skip_aligners
        
        script:
        """
        mkdir Hisat2Index
        hisat2_extract_exons.py $gtf > Hisat2Index/${gtf.baseName}.exon
        hisat2_extract_splice_sites.py $gtf > Hisat2Index/${gtf.baseName}.ss
        #hisat2_extract_splice_sites.py ${params.snp} > Hisat2Index/${gtf.baseName}.snp
        hisat2-build -p ${task.cpus} -f $fasta --exon Hisat2Index/${gtf.baseName}.exon --ss Hisat2Index/${gtf.baseName}.ss Hisat2Index/${fasta.baseName}
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no Hisat2 Index")
}

/*
 * PREPROCESSING - Build BWA index
 * NEED genome.fa
 */
if( params.bwa_index && !skip_bwa ){
    bwa_index = Channel
        .fromPath( params.bwa_index )
        .ifEmpty { exit 1, "bwa index not found: ${params.bwa_index}" }
}else if(params.fasta ){
    process MakeBWAIndex {
        label 'build_index'
        tag "bwa_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta

        output:
        file "BWAIndex/*" into bwa_index

        when:
        !skip_bwa && !skip_aligners
     
        script:
        """
        mkdir BWAIndex
        cd BWAIndex/
        bwa index -b ${task.cpus} -p ${fasta.baseName} -abwtsw ../$fasta
        cd ../
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no BWA Index")
}

/*
 * PREPROCESSING - Build STAR index
 * NEED genome.fa genes.gtf
 */
if( params.star_index && !skip_star){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}else if( params.fasta ){
    process MakeStarIndex {
        label 'build_index'
        tag "star_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        input:
        file fasta
        file gtf

        output:
        file "StarIndex" into star_index

        when:
        !skip_star && !skip_aligners 

        script:
        readLength = 50
        overhang = readLength - 1
        """
        mkdir StarIndex
        STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir StarIndex \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf \
        --sjdbOverhang $overhang 
        """
    }
}else {
   exit 1, println LikeletUtils.print_red("There is no STAR Index")
}
/*
========================================================================================
                                Step 1. QC------FastQC
========================================================================================
*/ 
process Fastp{
    tag "$sample_name"
    errorStrategy 'ignore'
    publishDir path: { params.skip_fastp ? params.outdir : "${params.outdir}/QC/fastp" },
             saveAs: { params.skip_fastp ? null : it }, mode: 'link'
        
    input:
    set sample_name, file(reads) from raw_data

    output:
    val sample_name into pair_id_fastqc, pair_id_tophat2, pair_id_hisat2, pair_id_bwa, pair_id_star 
    file "*_aligners.fastq" into fastqc_reads ,tophat2_reads, hisat2_reads, bwa_reads, star_reads
    file "*" into fastp_results

    when:
    !skip_aligners

    shell:
    skip_fastp = params.skip_fastp
    gzip = params.gzip
    if ( params.singleEnd ){
        if ( skip_fastp )  println LikeletUtils.print_purple("fastp is skipped")
        else println LikeletUtils.print_purple("fastp is on going for single-end data")
        filename = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners = sample_name + "_aligners.fastq"
        """
        if [ $gzip == "true" ]; then
            zcat ${reads} > ${sample_name}.fastq
        fi
        if [ $skip_fastp == "false" ]; then
            fastp -i ${sample_name}.fastq -o ${add_aligners} -j ${sample_name}_fastp.json -h ${sample_name}_fastp.html -w ${task.cpus}
        else
            mv ${sample_name}.fastq ${add_aligners}
        fi
        """
    } else {
        if ( skip_fastp )  println LikeletUtils.print_purple("fastp is skipped")
        else println LikeletUtils.print_purple("fastp is on going for pair-end data")
        filename = reads[0].toString() - ~/(_R[0-9])?(_[0-9])?(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners_1 = sample_name + "_1_aligners.fastq"
        add_aligners_2 = sample_name + "_2_aligners.fastq"
        """
        if [ $gzip == "true" ]; then
            zcat ${reads[0]} > ${sample_name}_1.fastq
            zcat ${reads[1]} > ${sample_name}_2.fastq
        fi
        if [ $skip_fastp == "false" ]; then  
            fastp -i ${sample_name}*1.fastq -o ${add_aligners_1} -I ${sample_name}*2.fastq -O ${add_aligners_2} -j ${sample_name}_fastp.json -h ${sample_name}_fastp.html -w ${task.cpus}
        else
            mv ${sample_name}*1.fastq ${add_aligners_1}
            mv ${sample_name}*2.fastq ${add_aligners_2}
        fi
        """
    } 
}
process Fastqc{
    tag "$sample_name"
    publishDir path: { params.skip_fastqc ? params.outdir : "${params.outdir}/QC" },
             saveAs: { params.skip_fastqc ? null : it }, mode: 'link'

    input:
    val sample_name from pair_id_fastqc
    file(reads) from fastqc_reads

    output:
    file "fastqc/*" into fastqc_results

    when:
    !skip_aligners && !params.skip_fastqc

    shell:
    skip_fastqc = params.skip_fastqc
    if ( params.singleEnd ){
        println LikeletUtils.print_purple("Fastqc is on going for single-end data")
        """
        mkdir fastqc
        fastqc -o fastqc --noextract ${reads}
        """       
    } else {
        println LikeletUtils.print_purple("Fastqc is on going for pair-end data")
        """
        mkdir fastqc   
        fastqc -o fastqc --noextract ${reads[0]}
        fastqc -o fastqc --noextract ${reads[1]}
        """      
    }
}
/*
========================================================================================
                            Step 2. Reads Mapping
========================================================================================
*/ 
process Tophat2Align {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/tophat2", mode: 'link', overwrite: true

    input:
    val sample_name from pair_id_tophat2
    file(reads) from tophat2_reads
    file index from tophat2_index.collect()
    file gtf

    output:
    file "*_tophat2.bam" into tophat2_bam
    file "*_log.txt" into tophat2_log
    
    when:
    !skip_tophat2 && !skip_aligners

    script:
    index_base = index[0].toString() - ~/(\.rev)?(\.\d)?(\.fa)?(\.bt2)?$/
    strand_info = params.stranded == "no" ? "fr-unstranded" : params.stranded == "reverse" ? "fr-secondstrand" : "fr-firststrand"
    if (params.singleEnd) {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Tophat2 in single-end mode")
        """
        tophat  -p ${task.cpus} \
                -G $gtf \
                -o $sample_name \
                --no-novel-juncs \
                --library-type $strand_info \
                $index_base \
                $reads &> ${sample_name}_log.txt
        mv $sample_name/accepted_hits.bam ${sample_name}_tophat2.bam
        """
    } else {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Tophat2 in paired-end mode")
        """
        tophat -p ${task.cpus} \
                -G $gtf \
                -o $sample_name \
                --no-novel-juncs \
                --library-type $strand_info \
                $index_base \
                ${reads[0]} ${reads[1]} &> ${sample_name}_log.txt
        mv $sample_name/accepted_hits.bam ${sample_name}_tophat2.bam
        """
    }
}

process Hisat2Align {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/hisat2", mode: 'link', overwrite: true

    input:
    val sample_name from pair_id_hisat2
    file(reads) from hisat2_reads
    file index from hisat2_index.collect()

    output:
    file "*_hisat2.bam" into hisat2_bam
    file "*_summary.txt" into hisat2_log

    when:
    !skip_hisat2 && !skip_aligners

    script:
    index_base = index[0].toString() - ~/(\.exon)?(\.\d)?(\.fa)?(\.gtf)?(\.ht2)?$/
    if (params.singleEnd) {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Hisat2 in single-end mode")
        """
        hisat2  -p ${task.cpus} --dta \
                -x $index_base \
                -U $reads \
                -S ${sample_name}_hisat2.sam &> ${sample_name}_hisat2_summary.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_hisat2.sam >${sample_name}_hisat2.bam
        rm *.sam
        """
    } else {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Hisat2 in paired-end mode")
        """
        hisat2  -p ${task.cpus} --dta \
                -x $index_base \
                -1 ${reads[0]} -2 ${reads[1]} \
                -S ${sample_name}_hisat2.sam &> ${sample_name}_hisat2_summary.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_hisat2.sam > ${sample_name}_hisat2.bam
        rm *.sam
        """
    }
}

process BWAAlign{
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/bwa", mode: 'link', overwrite: true
    
    input:
    val sample_name from pair_id_bwa
    file(reads) from bwa_reads
    file index from bwa_index.collect()

    output:
    file "*_bwa.bam" into bwa_bam
    file "*" into bwa_result

    when:
    !skip_bwa && !skip_aligners

    script:
    index_base = index[0].toString() - ~/(\.pac)?(\.bwt)?(\.ann)?(\.amb)?(\.sa)?(\.fa)?$/
    if (params.singleEnd) {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by BWA in single-end mode")
        """
        bwa aln -t ${task.cpus} \
                -f ${reads.baseName}.sai \
                $index_base \
                $reads
        bwa samse -f ${sample_name}_bwa.sam \
                $index_base \
                ${reads.baseName}.sai \
                $reads &> ${sample_name}_log.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_bwa.sam > ${sample_name}_bwa.bam
        rm *.sam
        """
    } else {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by BWA in paired-end mode")
        """
        bwa aln -t ${task.cpus} \
                -f ${reads[1].baseName}.sai \
                $index_base \
                ${reads[0]}
        bwa aln -t ${task.cpus} \
                -f ${reads[1].baseName}.sai \
                $index_base \
                ${reads[1]}
        bwa sampe -f ${sample_name}_bwa.sam \
                $index_base \
                ${reads[1].baseName}.sai ${reads[1].baseName}.sai \
                ${reads[0]} ${reads[1]} &> ${sample_name}_log.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_bwa.sam > ${sample_name}_bwa.bam
        rm *.sam
        """
    }
}

process StarAlign {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/star", mode: 'link', overwrite: true
    
    input:
    val sample_name from pair_id_star
    file(reads) from star_reads
    file star_index from star_index.collect()

    output:
    file "*_star.bam" into star_bam
    file "*.final.out" into star_log

    when:
    !skip_star && !skip_aligners

    script:
    if (params.singleEnd) {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by STAR in single-end mode")
        """
        STAR --runThreadN ${task.cpus} \
            --twopassMode Basic \
            --genomeDir $star_index \
            --readFilesIn $reads  \
            --outSAMtype BAM Unsorted \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFilterMultimapNmax 20 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix ${sample_name}  &> ${sample_name}_log.txt
        mv ${sample_name}Aligned.out.bam ${sample_name}_star.bam
        """
    } else {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by STAR in paired-end mode")
        """
        STAR --runThreadN ${task.cpus} \
            --twopassMode Basic \
            --genomeDir $star_index \
            --readFilesIn ${reads[0]} ${reads[1]}  \
            --outSAMtype BAM Unsorted \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFilterMultimapNmax 20 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix ${sample_name} &> ${sample_name}_log.txt
        mv ${sample_name}Aligned.out.bam ${sample_name}_star.bam
        """
    }
}
/*
========================================================================================
                        Step 3 Sort BAM file AND QC
========================================================================================
*/ 
if( params.aligners != "none"){
    Channel
        .from()
        .concat(tophat2_bam, hisat2_bam, bwa_bam, star_bam)
        .into {merge_bam_file; test_channel1}
}else{
    Channel
        .from()
        .concat(raw_bam)
        .into {merge_bam_file; test_channel1}
}

//test_channel1.subscribe{ println it }
/*
 * STEP 3-1 - Sort BAM file
*/
process Sort {
    tag "$sample_name"

    input:
    file bam_file from merge_bam_file

    output:
    file "*_sort*.{bam,bai}" into rename_bam_file
    file "*.bam" into bam_results

    script:
    sample_name = bam_file.toString() - ~/(\.bam)?$/
    output = sample_name + "_sort.bam"
    keep_unique = (params.mapq_cutoff).toInteger() 
    if (!params.skip_sort){
        println LikeletUtils.print_purple("Samtools sorting the bam files now")
        if(keep_unique) println LikeletUtils.print_purple("Keeping unique reads")
        """
        if [ "$keep_unique" -gt "0" ]; then
            samtools view -bq $keep_unique $bam_file | samtools sort -@ ${task.cpus} -O BAM -o $output -
        else
            samtools sort -@ ${task.cpus} -O BAM -o $output $bam_file
        fi
        samtools index -@ ${task.cpus} $output
        """
    } else {
        println LikeletUtils.print_purple("The step of samtools sort is skipped")
        if(keep_unique) println LikeletUtils.print_purple("Keeping unique reads")
        """
        for bam_file in *.bam
        do
        {
            if [ "$keep_unique" -gt "0" ]; then
                samtools view -bq $keep_unique $bam_file > $output
            else
                mv $bam_file $output
            fi
            samtools index -@ ${task.cpus} $output
        }
        done
        """
    }
}

process RenameByDesignfile{
    publishDir "${params.outdir}/alignment/samtoolsSort/", mode: 'link', overwrite: true

    input:
    file (reads) from rename_bam_file.collect()
    file designfile  // designfile:filename,control_treated,input_ip

    output:
    file "*.{input,ip}_*.{bam,bai}" into sort_bam
    file ("formatted_designfile.txt") into formatted_designfile
    file "*" into rename_results

    script:
    println LikeletUtils.print_purple("Rename the files for downstream analysis")
    aligners_name = params.aligners
    """
    # Windows and linux newline ^M conversion
    cat $designfile | dos2unix |sed '1s/.*/Sample_ID,input_FileName,ip_FileName,Group/g' > formatted_designfile.txt
    # Check the consistency of designfile and comparefile
    if [ "$comparefile" != "false" ]; then 
        ## get groups' name in comparefile
        cat $comparefile | dos2unix | awk -F "_vs_" '{print \$1"\\n"\$2}' | sort | uniq > tmp.compare.group
        ## get groups' name in designfile
        awk -F, 'NR>1{print \$4}' formatted_designfile.txt | sort | uniq > tmp.design.group
        intersection_num=\$(join tmp.compare.group tmp.design.group | wc -l)
        if [[ \$intersection_num  != \$(cat tmp.compare.group| wc -l) ]] ;then 
            echo "The groups' name of comparefile and designfile are inconsistent."
            echo "Please check your comparefile: "$comparefile
            echo "The groups' name of comparefile in designfile: "\$(join tmp.compare.group tmp.design.group)
            exit 1
        fi
        rm tmp.compare.group tmp.design.group
    fi
    bash $baseDir/bin/rename.sh $aligners_name formatted_designfile.txt
    bash $baseDir/bin/rename_for_resume.sh formatted_designfile.txt
    """
}
/*
 * STEP 3-2 - RSeQC analysis
*/
process RSeQC {
    publishDir "${params.outdir}/QC/rseqc" , mode: 'copy', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }    
    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc from sort_bam.collect()
    file bed12 from bed_rseqc.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results
    file "*.bam_stat.txt" into bam_stat_for_normlization
    script:

    """    
    bash $baseDir/bin/rseqc.sh $bed12 ${task.cpus}
    """
}

process CreateBigWig {
    publishDir "${params.outdir}/QC/rseqc/", mode: 'link', overwrite: true ,
        saveAs: {filename ->
            if (filename.indexOf("bigwig") > 0) "bigwig/$filename"
        }

    input:
    file bam from sort_bam.collect()

    output:
    file "*.bigwig" into bigwig_for_genebody

    when:
    !params.skip_createbigwig

    script:
    """
    bash $baseDir/bin/create_bigwig.sh ${task.cpus}
    """
}

process GenebodyCoverage {
    publishDir "${params.outdir}/QC/rseqc" , mode: 'link', overwrite: true, 
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_rseqc && !params.skip_genebody_coverage

    input:
    file bigwig from bigwig_for_genebody
    file bed12 from bed_genebody_coverage.collect()

    output:
    file "*.{txt,pdf,r}" into genebody_coverage_results

    script:
    """
    bash $baseDir/bin/geneBody_coverage2.sh $bed12 ${task.cpus}
    """
}
Channel
    .from()
    .concat( tophat2_log, hisat2_log, star_log, fastp_results, fastqc_results, rseqc_results)
    .into{ arranged_qc; qc_results_for_report }

process multiqc{
    label 'report'
    publishDir "${params.outdir}/Report/QCReadsReport" , mode: 'link', overwrite: true
    
    when:
    !params.skip_qc

    input:
    file arranged_qc from arranged_qc.collect()

    output:
    file "multiqc*" into multiqc_results

    script:
    aligner = params.aligners
    """
    multiqc -n multiqc_$aligner .
    """
}
/*
========================================================================================
                            Step 4 Peak Calling
========================================================================================
*/ 
/*
 * STEP 4 - 1  Peak Calling------MetPeak, MACS2, MATK
*/
process Metpeak {
    label 'peak_calling'
    publishDir "${params.outdir}/peakCalling/metpeak", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file gtf
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*" into metpeak_results
    file "metpeak*_normalized.bed" into metpeak_nomarlized_bed

    when:
    !params.skip_metpeak && !params.skip_peakCalling

    script: 
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    if( flag_peakCallingbygroup ){
        println LikeletUtils.print_purple("Peak Calling performed by MeTPeak in group mode")
    }else{
        println LikeletUtils.print_purple("Peak Calling performed by MeTPeak in independent mode")
    }
    """
    Rscript $baseDir/bin/MeTPeak.R $formatted_designfile $gtf ${task.cpus} $flag_peakCallingbygroup;
    """
}

process Macs2{
    label 'peak_calling'
    publishDir "${params.outdir}/peakCalling/macs2", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "macs2*.{xls,narrowPeak,summits}" into macs2_results
    file "macs2*_normalized.bed" into macs2_nomarlized_bed

    when:
    !params.skip_macs2 && !params.skip_peakCalling

    script:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    if( flag_peakCallingbygroup ){
        println LikeletUtils.print_purple("Peak Calling performed by Macs2 in group mode")
    }else{
        println LikeletUtils.print_purple("Peak Calling performed by Macs2 in independent mode")
    }
    """
    bash $baseDir/bin/macs2.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup;
    """ 
}

process MATKpeakCalling {
    label 'peak_calling'
    publishDir "${params.outdir}/peakCalling/MATK", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file gtf
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*" into matk_results
    file "MATK*_normalized.bed" into matk_nomarlized_bed

    when:
    !params.skip_matk && !params.skip_peakCalling

    script:
    matk_jar = params.matk_jar
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    if( flag_peakCallingbygroup ){
        println LikeletUtils.print_purple("Peak Calling performed by MATK in group mode")
    }else{
        println LikeletUtils.print_purple("Peak Calling performed by MATK in independent mode")
    }
    """
    export OMP_NUM_THREADS=${task.cpus}
    bash $baseDir/bin/MATK_peakCalling.sh $matk_jar $formatted_designfile $flag_peakCallingbygroup
    """    
}

process Meyer{
    label 'peak_calling'
    publishDir "${params.outdir}/peakCalling/meyer", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file chromsizesfile from chromsizesfile

    output:
    file "meyer*.bed" into meyer_results
    file "meyer*_normalized.bed" into meyer_nomarlized_bed

    when:
    !params.skip_meyer && !params.skip_peakCalling

    shell:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    if( flag_peakCallingbygroup ){
        println LikeletUtils.print_purple("Peak Calling performed by Meyer in group mode")
    }else{
        println LikeletUtils.print_purple("Peak Calling performed by Meyer in independent mode")
    }
    '''
    cp !{baseDir}/bin/meyer.py ./
    awk '{print $1}' !{chromsizesfile} > chrName.txt
    mkdir genomebin
    bedtools makewindows -g !{chromsizesfile} -w 25 > genome.bin25.bed
    peak_windows_number=$(wc -l genome.bin25.bed| cut -d " " -f 1)
    awk '{print "cat genome.bin25.bed | grep "$1" > genomebin/"$1".bin25.bed"}' chrName.txt | xargs -iCMD -P!{task.cpus} bash -c CMD
    bash !{baseDir}/bin/meyer_peakCalling.sh !{formatted_designfile} chrName.txt genomebin/ ${peak_windows_number} !{task.cpus} !{flag_peakCallingbygroup}
    ''' 
}
/*
========================================================================================
                        Step 5 Differential expression analysis
========================================================================================
*/
process HtseqCount{
    label 'analysis'
    publishDir "${params.outdir}/expressionAnalysis/htseq-count", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file gtf

    output:
    file "*input*.count" into htseq_count_input_to_deseq2, htseq_count_input_to_edgeR, htseq_count_input_to_arrange
    file "expression.matrix" into htseq_results

    when:
    !params.skip_expression

    script:
    println LikeletUtils.print_purple("Generate gene expression matrix by htseq-count and Rscript")
    strand_info = params.stranded == "no" ? "no" : params.stranded == "reverse" ? "reverse" : "yes"
    """
    bash $baseDir/bin/htseq_count.sh $gtf $strand_info ${task.cpus} 
    Rscript $baseDir/bin/get_htseq_matrix.R $formatted_designfile  ${task.cpus} 
    """ 
}

process DESeq2{
    label 'analysis'
    tag "$compare_str"

    publishDir "${params.outdir}/expressionAnalysis/DESeq2", mode: 'link', overwrite: true

    input:
    file reads_count_input from htseq_count_input_to_deseq2
    file formatted_designfile from formatted_designfile.collect()
    val compare_str from compareLines_for_DESeq2

    output:
    file "DESeq2*.csv" into deseq2_results
    
    when:
    !params.skip_deseq2 && !params.skip_expression
    
    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by DESeq2 ($compare_str)")
    """
    Rscript $baseDir/bin/DESeq2.R $formatted_designfile $compare_str
    """ 
}

process EdgeR{
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/expressionAnalysis/edgeR", mode: 'link', overwrite: true

    input:
    file reads_count_input from htseq_count_input_to_edgeR
    file formatted_designfile from formatted_designfile.collect()
    val compare_str from compareLines_for_edgeR

    output:
    file "edgeR*.csv" into edgeR_results
    
    when:
    !params.skip_edger && !params.skip_expression

    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by EdgeR ($compare_str)")
    """
    Rscript $baseDir/bin/edgeR.R $formatted_designfile $compare_str
    """ 
}

process Cufflinks{
    label 'analysis'
    publishDir "${params.outdir}/expressionAnalysis/cufflinks", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file gtf
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "cuffdiff*" into cufflinks_results

    when:
    !params.skip_cufflinks && !params.skip_expression

    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by Cufflinks")
    """
    bash $baseDir/bin/cufflinks.sh $formatted_designfile $gtf ${task.cpus}
    """
}
/*
========================================================================================
                        Step 6 Merge Peak AND Analysis
========================================================================================
*/
/*
 * STEP 5-1 Merge Peak
*/

Channel
    .from()
    .concat(metpeak_nomarlized_bed, macs2_nomarlized_bed, matk_nomarlized_bed, meyer_nomarlized_bed)
    .into {merged_bed ; bed_for_annotation}

process PeakMerge {
    label 'analysis'
    publishDir "${params.outdir}/peakCalling/mergedBed", mode: 'link', overwrite: true,
        saveAs: {filename ->
            if (filename.indexOf("bed") > 0) "$filename"
        }
    
    input:
    file peak_bed from merged_bed.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*merged*.bed" into merge_result
    file "*merged_group*.bed" into group_merged_bed
    file "*_merged_allpeaks.bed" into all_merged_bed

    script:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    peakCalling_tools_count = (params.skip_metpeak ? 0 : 1).toInteger() + (params.skip_macs2 ? 0 : 1).toInteger() + (params.skip_matk ? 0 : 1).toInteger() + (params.skip_meyer ? 0 : 1).toInteger()
    peakMerged_mode = params.peakMerged_mode
    if ( peakMerged_mode == "rank" )  
        println LikeletUtils.print_purple("Start merge peaks by RobustRankAggreg")
    else if ( peakMerged_mode == "mspc" )  
        println LikeletUtils.print_purple("Start merge peaks by MSPC")
    else
        println LikeletUtils.print_purple("Start merge peaks by " + peakMerged_mode )
    """
    cp ${baseDir}/bin/normalize_peaks.py ./
    if [ ${peakMerged_mode} == "rank" ]; then 
        cp $baseDir/bin/merge_peaks_by_rank.R ./
        bash $baseDir/bin/merge_peaks_by_rank.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count
    elif [ ${peakMerged_mode} == "mspc" ]; then
        bash $baseDir/bin/merge_peaks_by_mspc.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count mspc_results
    elif [ ${peakMerged_mode} == "macs2" ]||[ ${peakMerged_mode} == "metpeak" ]||[ ${peakMerged_mode} == "MATK" ]||[ ${peakMerged_mode} == "meyer" ]; then 
        bash $baseDir/bin/merge_peaks_by_bedtools.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count $peakMerged_mode
    else
        echo -e "Please check your value of peakMerged_mode: $peakMerged_mode"
    fi
    """
}

Channel
    .from()
    .concat( group_merged_bed, all_merged_bed )
    .into { annotate_collection; motif_collection; bed_collect_for_arrange_results}

process BedAnnotated{
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/AnnotatedPeaks", mode: 'link', overwrite: true
    
    input:
    file all_bed from bed_for_annotation.collect()
    file annotate_file from annotate_collection.collect()
    file formatted_designfile from formatted_designfile.collect()
    file fasta
    file gtf

    output:
    file "annotatedby{xy,homer}/*" into annotation_results,annotation_results_2
    file "annotatedbyxy/*merged_allpeaks.anno.txt" into methylation_annotaion_file
    
    when:
    !params.skip_annotation

    script:
    annotated_script_dir = baseDir + "/bin"
    //Skip Peak Calling Tools Setting
    """
    # Annotation Peaks
    cp ${annotated_script_dir}/intersec.pl ./
    cp ${annotated_script_dir}/m6A_annotate_forGTF_xingyang2.pl ./
    bash ${baseDir}/bin/annotation.sh ${fasta} ${gtf} ${task.cpus}  
    """
}
process MotifSearching {
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/motif", mode: 'link', overwrite: true
    
    input:
    file group_bed from motif_collection.collect()
    file formatted_designfile from formatted_designfile.collect()
    file fasta
    file gtf

    output:
    file "*_{dreme,homer}" into motif_results,motif_results_2

    when:
    !params.skip_motif

    script:
    motif_file_dir = baseDir + "/bin"
    println LikeletUtils.print_purple("Motif analysis is going on by DREME and Homer")
    """
    cp ${motif_file_dir}/m6A_motif.meme ./
    bash $baseDir/bin/motif_searching.sh $fasta $gtf m6A_motif.meme ${task.cpus}
    """
}

process QCPeaksReport {
    label 'report'
    publishDir "${params.outdir}/Report/QCPeaksReport", mode: 'link', overwrite: true
    
    input:
    file motif from motif_results.collect()
    file annotation_files from annotation_results.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*.{html,pdf}" into qcPeaksReport

    script:
    peakMerged_mode = params.peakMerged_mode 
    peakCalling_mode = params.peakCalling_mode
    qcPeaksRData = "QCPeakPlot.RData"
    """
    cp $baseDir/bin/QC_Peaks_Report.rmd ./
    Rscript $baseDir/bin/QC_Peaks_Report.R $formatted_designfile $peakMerged_mode $peakCalling_mode $qcPeaksRData
    R -e "load(\\"$qcPeaksRData\\");rmarkdown::render('QC_Peaks_Report.rmd',output_file='QC_Peaks_Report_${peakMerged_mode}.html')"
    """
}

process PeaksQuantification{
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/m6AQuantification", mode: 'link', overwrite: true
    
    input:
    file merged_bed from all_merged_bed.collect()
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file annotation_file from methylation_annotaion_file.collect()
    file gtf

    output:
    file "*.{matrix,count}" into quantification_results, quantification_matrix

    when:
    !params.skip_peakCalling

    script:
    matk_jar = params.matk_jar
    methylation_analysis_mode = params.methylation_analysis_mode
    if ( methylation_analysis_mode == "Wilcox-test" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by bedtools")
    else if ( methylation_analysis_mode == "QNB" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by QNB")
    else if ( methylation_analysis_mode == "MATK" )
        println LikeletUtils.print_purple("Generate m6A quantification matrix by MATK")
    else if ( methylation_analysis_mode == "edgeR" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by edgeR")
    else if ( methylation_analysis_mode == "DESeq2" )
        println LikeletUtils.print_purple("Generate m6A quantification matrix by DESeq2")
    """
    if [ ${methylation_analysis_mode} == "Wilcox-test" ]; then 
        # PeaksQuantification by bedtools
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/bedtools_quantification.R $formatted_designfile bam_stat_summary.txt
    elif [ ${methylation_analysis_mode} == "QNB" ]||[ ${methylation_analysis_mode} == "MeTDiff" ]; then 
        # PeaksQuantification by QNB
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/QNB_quantification.R $formatted_designfile ${task.cpus}
    elif [ ${methylation_analysis_mode} == "MATK" ]; then 
        # PeaksQuantification by MATK
        export OMP_NUM_THREADS=${task.cpus}
        bash $baseDir/bin/MATK_quantification.sh $matk_jar $gtf $formatted_designfile ${merged_bed} 1
    elif [ ${methylation_analysis_mode} == "DESeq2" ]||[ ${methylation_analysis_mode} == "edgeR" ]; then 
        # PeaksQuantification by LRT
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/bedtools_quantification.R $formatted_designfile bam_stat_summary.txt
    fi
    head -1 *_quantification.matrix |sed 's/^\\t//'  |awk -F "\\t" '{print "ID\\tGene_symbol\\t"\$0}' > tmp.header.file
    sed '1d' *_quantification.matrix | sort > tmp.quantification.file
    awk 'BEGIN{FS="\\t";OFS="\\t"}{print \$4,\$15,\$11}' ${annotation_file} | sort | join -t \$'\t' -e 'NA' -a1 -o 1.1 -o 2.2 -o 2.3 tmp.quantification.file - >  tmp.annotation.file
    join -t \$'\t' tmp.annotation.file tmp.quantification.file | cat tmp.header.file - > *_quantification.matrix 
    """
}

process diffm6APeak{
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/m6AAnalysis/diffm6A", mode: 'link', overwrite: true
    
    input:
    //file peak_bed from group_merged_bed.collect()
    file merged_bed from all_merged_bed.collect()
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file count_matrix from quantification_matrix.collect()
    file gtf
    val compare_str from compareLines_for_diffm6A

    output:
    file "*diffm6A*.txt" into diffm6A_results

    when:
    !params.skip_diffpeakCalling && params.comparefile

    script:
    matk_jar = params.matk_jar
    methylation_analysis_mode = params.methylation_analysis_mode
    if ( methylation_analysis_mode == "Wilcox-test" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by bedtools")
    else if ( methylation_analysis_mode == "QNB" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by QNB")
    else if ( methylation_analysis_mode == "MATK" )
        println LikeletUtils.print_purple("Differential m6A analysis is going on by MATK")
    else if ( methylation_analysis_mode == "edgeR" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by edgeR")
    else if ( methylation_analysis_mode == "DESeq2" )
        println LikeletUtils.print_purple("Differential m6A analysis is going on by DESeq2")
    """
    if [ ${methylation_analysis_mode} == "Wilcox-test" ]; then 
        Rscript $baseDir/bin/bedtools_diffm6A.R $formatted_designfile bedtools_quantification.matrix $compare_str
    elif [ ${methylation_analysis_mode} == "QNB" ]; then 
        Rscript $baseDir/bin/QNB_diffm6A.R $formatted_designfile ${merged_bed} $compare_str     
    elif [ ${methylation_analysis_mode} == "MeTDiff" ]; then
        Rscript $baseDir/bin/MeTDiff_diffm6A.R $formatted_designfile $compare_str 
    elif [ ${methylation_analysis_mode} == "MATK" ]; then 
        export OMP_NUM_THREADS=${task.cpus}
        bash $baseDir/bin/MATK_diffm6A.sh $matk_jar $formatted_designfile $gtf $compare_str $merged_bed
    elif [ ${methylation_analysis_mode} == "edgeR" ]; then
        Rscript $baseDir/bin/GLM_edgeR_DM.R $formatted_designfile $compare_str 
    elif [ ${methylation_analysis_mode} == "DESeq2" ]; then 
        Rscript $baseDir/bin/GLM_DESeq2_DM.R $formatted_designfile $compare_str ${task.cpus}
    fi

    """ 
}

process SingleNucleotidePrediction{
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/m6APredictionSites", mode: 'link', overwrite: true
    
    input:
    file peak_bed from group_merged_bed.collect()
    file group_bed from all_merged_bed.collect()
    file formatted_designfile from formatted_designfile.collect()
    file bam_bai_file from sort_bam.collect()
    file fasta
    file gtf

    output:
    file "m6A_sites*.bed" into prediction_results

    when:
    !params.skip_m6Aprediction

    script:
    matk_jar = params.matk_jar
    println LikeletUtils.print_purple("SignleNucleotide Prediction analysis is going on by MATK")
    """
    export OMP_NUM_THREADS=${task.cpus}
    bash $baseDir/bin/m6Aprediction.sh $matk_jar $formatted_designfile $fasta $gtf
    """
}


Channel
    .from()
    .concat( quantification_results, motif_results_2, diffm6A_results, 
        htseq_count_input_to_arrange, 
        annotation_results_2, prediction_results, bed_collect_for_arrange_results,
        multiqc_results, deseq2_results, edgeR_results, cufflinks_results
    )
    .set{ results_arrange }

process DiffReport {
    label 'report'
    publishDir "${params.outdir}/Report" , mode: 'link', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf(".html") > 0)  "diffReport/$filename"
                 else if (filename.indexOf(".pdf") > 0)  "diffReport/$filename"
                 else "ReportRData/$filename"
        }        
    input:
    file results from results_arrange.collect()
    file formatted_designfile from formatted_designfile.collect()
    val compare_info from compareLines_for_arranged_result.collect()
    
    output:
    file "*.m6APipe" into m6APipe_result
    file "*.{html,pdf}" into diffReport_result

    // when:
    // !params.skip_annotation && !params.skip_expression && !skip_diffpeakCalling

    script:
    methylation_analysis_mode = params.methylation_analysis_mode
    expression_analysis_mode = params.expression_analysis_mode
    peakMerged_mode = params.peakMerged_mode
    diffReportRData = "DiffReport.RData"
    """
    cp $baseDir/bin/DiffReport.rmd ./
    if [ "$compare_info" != "[two_group]" ]; then
        echo $compare_info | sed 's/^\\[//g' | sed 's/\\]\$//g' | sed s/[[:space:]]//g > compare_info
    else
        echo \$(awk 'BEGIN{FS=","}NR>1{print \$4}' $formatted_designfile |sort|uniq|awk 'NR==1{printf \$0"_vs_"}NR==2{print \$0}') > compare_info
    fi
    Rscript $baseDir/bin/arranged_results.R $formatted_designfile compare_info $methylation_analysis_mode $expression_analysis_mode $peakMerged_mode
    Rscript $baseDir/bin/DiffReport.R *.m6APipe $diffReportRData
    R -e "load(\\"$diffReportRData\\");rmarkdown::render('DiffReport.rmd',output_file='DiffReport_${peakMerged_mode}_${methylation_analysis_mode}_${expression_analysis_mode}.html')"
    rm Rplots.pdf
    """
}

process CreateIGVjs {
    label 'report'
    publishDir "${params.outdir}/Report" , mode: 'link', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf(".html") > 0)  "Igv_js/$filename"
                 else if (filename.indexOf(".pdf") > 0)  "Igv_js/$filename"
                 else "Igv_js/$filename"
        }        
    input:
    file m6APipe_result from m6APipe_result
    file fasta 
    file gtf
    file formatted_designfile from formatted_designfile.collect()
    file group_bed from group_merged_bed.collect()
    file all_bed from all_merged_bed.collect()
    file bigwig from bigwig_for_genebody.collect()
    
    output:
    file "*" into igv_js

    script:    
    igv_fasta = fasta.baseName.toString() + ".igv.fa"
    igv_gtf = gtf.baseName.toString() + ".igv.gtf"
    merged_allpeaks_igvfile = all_bed.baseName.toString() + ".igv.bed"
    """
    ls -l $fasta | awk -F "> " '{print "ln "\$2" ./'$igv_fasta'"}' | bash
    ls -l $gtf | awk -F "> " '{print "ln "\$2" ./'$igv_gtf'"}' | bash
    ls -l $m6APipe_result | awk '{print "ln "\$11" initial.m6APipe"}' | bash
    ls -l $group_bed $all_bed | awk '{sub(".bed\$",".igv.bed",\$9);print "ln "\$11,\$9}' | bash
    ls -l $bigwig | awk '{sub(".bigwig\$",".igv.bigwig",\$9);print "ln "\$11,\$9}' | bash
    samtools faidx $igv_fasta
    bash $baseDir/bin/create_IGV_js.sh $fasta $gtf $merged_allpeaks_igvfile $formatted_designfile
    """
}

/*
Working completed message
 */
workflow.onComplete {
    println LikeletUtils.print_green("=================================================")
    println LikeletUtils.print_green("Cheers! m6APipe from SYSUCC run Complete!")
    println LikeletUtils.print_green("=================================================")
    //email information
    if (params.mail) {
        recipient = params.mail
        def subject = 'My m6Aseq-SYSUCC execution'

        def  msg = """\
            RNAseq-SYSUCC execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """.stripIndent()

        sendMail(to: recipient, subject: subject, body: msg)
    }
}
workflow.onError {
   println LikeletUtils.print_yellow("Oops... Pipeline execution stopped with the following message: ")+LikeletUtils.print_red(workflow.errorMessage)
}
