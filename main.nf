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
unstranded = params.unstranded ? true : false


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
if( params.chromsizesfile ) {
    chromsizesfile = file(params.chromsizesfile)
    if( !chromsizesfile.exists() ) exit 1, LikeletUtils.print_red("Chromsizes file not found: ${params.chromsizesfile}")
}else{
    chromsizesfile = ""
}
if( params.comparefile == "two_group" ){
    compareLines = Channel.from("two_group")
}else if( params.comparefile){
    File comparefile = new File(params.comparefile)
    if( !comparefile.exists() ) exit 1, print_red("Compare file not found: ${params.comparefile}")
    compareLines = Channel.from(comparefile.readLines())
}else {
    compareLines = Channel.from("")
}
compareLines.into{compareLines_for_Deseq2; compareLines_for_edgeR; compareLines_for_diffm6A; compareLines_for_arranged_result }

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
}
else{
    exit 1, LikeletUtils.print_red("Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'tophat2', 'bwa'")
}

/*
 * Create a channel for input read files
 */
if( params.readPaths && !skip_aligners ){
    if( params.singleEnd ){
        Channel
            .fromFilePairs( "${params.readPaths}/*.fastq", size: 1 ) 
            .ifEmpty { exit 1, LikeletUtils.print_red("readPaths was empty - no fastq files supplied: ${params.readPaths}")}
            .into{ raw_data; raw_bam }
    }
    else if ( !params.singleEnd ){
        Channel
            .fromFilePairs( "${params.readPaths}/*_{1,2}.fastq", size: 2 ) 
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
log.info LikeletUtils.sysucc_ascii()
log.info LikeletUtils.print_purple("============You are running m6APipe with the following parameters===============")
log.info LikeletUtils.print_purple("Checking parameters ...")

log.info LikeletUtils.print_yellow("=====================================Reads types================================")
log.info (LikeletUtils.print_yellow("SingleEnd :                     ") + LikeletUtils.print_green(params.singleEnd))
log.info (LikeletUtils.print_yellow("Strand :                        ") + LikeletUtils.print_green(params.strand))

log.info LikeletUtils.print_yellow("====================================Mode selected==============================")
log.info (LikeletUtils.print_yellow("aligners :                      ") + LikeletUtils.print_green(params.aligners))
log.info (LikeletUtils.print_yellow("peakCalling_mode :              ") + LikeletUtils.print_green(params.peakCalling_mode))
log.info (LikeletUtils.print_yellow("mergepeak_mode :                ") + LikeletUtils.print_green(params.mergepeak_mode))
log.info (LikeletUtils.print_yellow("quantification_mode :           ") + LikeletUtils.print_green(params.quantification_mode))

log.info LikeletUtils.print_yellow("==================================Input files selected==========================")
log.info (LikeletUtils.print_yellow("Reads Path:                     ") + LikeletUtils.print_green(params.readPaths))
log.info (LikeletUtils.print_yellow("fasta file :                    ") + LikeletUtils.print_green(params.fasta))
log.info (LikeletUtils.print_yellow("Gtf file :                      ") + LikeletUtils.print_green(params.gtf))
log.info (LikeletUtils.print_yellow("Bed12 file :                    ") + LikeletUtils.print_green(params.bed12))
log.info (LikeletUtils.print_yellow("Design file :                   ") + LikeletUtils.print_green(params.designfile))
log.info (LikeletUtils.print_yellow("Compare file :                  ") + LikeletUtils.print_green(params.comparefile))

log.info LikeletUtils.print_yellow("==================================Skip model selected==========================")
log.info (LikeletUtils.print_yellow("Skip samtools sort :            ") + LikeletUtils.print_green(params.skip_sort))
log.info (LikeletUtils.print_yellow("Skip expression analysis :      ") + LikeletUtils.print_green(params.skip_expression))
log.info (LikeletUtils.print_yellow("Skip peakCalling :              ") + LikeletUtils.print_green(params.skip_peakCalling))
log.info (LikeletUtils.print_yellow("Skip diffpeakCalling :          ") + LikeletUtils.print_green(params.skip_diffpeakCalling))
log.info (LikeletUtils.print_yellow("Skip annotation :               ") + LikeletUtils.print_green(params.skip_annotation))
log.info (LikeletUtils.print_yellow("Skip qc :                       ") + LikeletUtils.print_green(params.skip_qc))

log.info LikeletUtils.print_yellow("==================================Output files directory========================")
log.info (LikeletUtils.print_yellow("Output directory :              ") + LikeletUtils.print_green(params.outdir))


/*
========================================================================================
                             check or build the index
========================================================================================
*/ 
/*
 * PREPROCESSING - Build BED12 file
 * NEED gtf.file
 */
if( params.bed12 ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, LikeletUtils.print_red("BED12 annotation file not found: ${params.bed12}") }
        .into {bed_rseqc; bed_genebody_coverage}
}else if( params.gtf && !params.bed12 ){
    process makeBED12 {
        tag "gtf2bed12"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf
        
        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script:      
        """
        bash ${baseDir}/bin/gtf2bed12.sh $gtf
        """        
    }
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
        hisat2-build -p ${task.cpus} -f $fasta --snp Hisat2Index/${gtf.baseName}.snp --exon Hisat2Index/${gtf.baseName}.exon --ss Hisat2Index/${gtf.baseName}.ss Hisat2Index/${fasta.baseName}
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
    publishDir path: { params.skip_fastp ? params.outdir : "${params.outdir}/QC/fastp" },
             saveAs: { params.skip_fastp ? null : it }, mode: 'link'
        
    input:
    set sample_name, file(reads) from raw_data

    output:
    val sample_name into pair_id_fastqc, pair_id_tophat2, pair_id_hisat2, pair_id_bwa, pair_id_star 
    file "*.fastq" into fastqc_reads ,tophat2_reads, hisat2_reads, bwa_reads, star_reads
    file "*" into fastp_results

    when:
    !skip_aligners

    shell:
    skip_fastp = params.skip_fastp
    if ( params.singleEnd ){
        if ( skip_fastp )  println LikeletUtils.print_purple("fastp is skipped")
        else println LikeletUtils.print_purple("fastp is on going for single-end data")
        filename = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners = reads.toString() - ~/(\.fq)?(\.fastq)?$/ + "_aligners.fastq"
        """
        if [ $skip_fastp == "false" ]; then
            fastp -i ${reads} -o ${add_aligners} 
        else
            mv ${reads} ${add_aligners}
        fi
        """
    } else {
        if ( skip_fastp )  println LikeletUtils.print_purple("fastp is skipped")
        else println LikeletUtils.print_purple("fastp is on going for pair-end data")
        filename = reads[0].toString() - ~/(_R[0-9])(_[0-9])?(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners_1 = reads[0].toString() - ~/(\.fq)?(\.fastq)?$/ + "_aligners.fastq"
        add_aligners_2 = reads[1].toString() - ~/(\.fq)?(\.fastq)?$/ + "_aligners.fastq"
        """
        if [ $skip_fastp == "false" ]; then  
            fastp -i ${reads[0]} -o ${add_aligners_1} -I ${reads[1]} -O ${add_aligners_2}
        else
            mv ${reads[0]} ${add_aligners_1}
            mv ${reads[1]} ${add_aligners_2}
        fi
        """
    }
}
process Fastqc{
    tag "$sample_name"
    publishDir path: { params.skip_fastqc ? params.outdir : "${params.outdir}/QC/fastqc" },
             saveAs: { params.skip_fastqc ? null : it }, mode: 'link'

    input:
    val sample_name from pair_id_fastqc
    file(reads) from fastqc_reads

    output:
    file "*" into fastqc_results

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
    tag "$sample_name"
    publishDir "${params.outdir}/aligners/tophat2", mode: 'link', overwrite: true

    input:
    val sample_name from pair_id_tophat2
    file(reads) from tophat2_reads
    file index from tophat2_index.collect()
    file gtf

    output:
    file "*_tophat2.bam" into tophat2_bam
    file "*" into tophat2_result
    
    when:
    !skip_tophat2 && !skip_aligners

    script:
    index_base = index[0].toString() - ~/(\.rev)?(\.\d)?(\.fa)?(\.bt2)?$/
    strand_str = unstranded ? "fr-unstranded" : "fr-firststrand"
    if (params.singleEnd) {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Tophat2 in single-end mode")
        """
        tophat  -p ${task.cpus} \
                -G $gtf \
                -o $sample_name \
                --no-novel-juncs \
                --library-type $strand_str \
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
                --library-type $strand_str \
                $index_base \
                ${reads[0]} ${reads[1]} &> ${sample_name}_log.txt
        mv $sample_name/accepted_hits.bam ${sample_name}_tophat2.bam
        """
    }
}

process Hisat2Align {
    tag "$sample_name"
    publishDir "${params.outdir}/aligners/hisat2", mode: 'link', overwrite: true

    input:
    val sample_name from pair_id_hisat2
    file(reads) from hisat2_reads
    file index from hisat2_index.collect()

    output:
    file "*_hisat2.bam" into hisat2_bam
    file "*" into hisat2_result

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
                -S ${sample_name}_hisat2.sam 2> ${sample_name}_hisat2_summary.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_hisat2.sam >${sample_name}_hisat2.bam
        rm *.sam
        """
    } else {
        println LikeletUtils.print_purple("Initial reads mapping of " + sample_name + " performed by Hisat2 in paired-end mode")
        """
        hisat2  -p ${task.cpus} --dta \
                -x $index_base \
                -1 ${reads[0]} -2 ${reads[1]} \
                -S ${sample_name}_hisat2.sam 2> ${sample_name}_hisat2_summary.txt
        samtools view -@ ${task.cpus} -h -bS ${sample_name}_hisat2.sam > ${sample_name}_hisat2.bam
        rm *.sam
        """
    }
}

process BWAAlign{
    tag "$sample_name"
    publishDir "${params.outdir}/aligners/bwa", mode: 'link', overwrite: true
    
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
    tag "$sample_name"
    publishDir "${params.outdir}/aligners/star", mode: 'link', overwrite: true
    
    input:
    val sample_name from pair_id_star
    file(reads) from star_reads
    file star_index from star_index.collect()

    output:
    file "*_star.bam" into star_bam
    file "*.final.out" into star_result

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
Channel
    .from()
    .concat(tophat2_bam, hisat2_bam, bwa_bam, star_bam, raw_bam)
    .into {merge_bam_file; test_channel1}
//test_channel1.subscribe{ println it }
/*
 * STEP 3-1 - Sort BAM file
*/
process RenameByDesignfile{
    //publishDir "${params.outdir}/Sample_rename", mode: 'link', overwrite: true

    input:
    file (reads) from merge_bam_file.collect()
    file designfile  // designfile:filename,control_treated,input_ip

    output:
    file ("*.bam") into rename_bam_file
    file ("formatted_designfile.txt") into formatted_designfile
    
    when:
    true

    script:
    println LikeletUtils.print_purple("Rename the files for downstream analysis")
    aligners_name = params.aligners
    """
    #Windows and linux newline ^M conversion
    cat $designfile | dos2unix |sed '1s/.*/Sample_ID,input_FileName,ip_FileName,Group/g' > formatted_designfile.txt
    bash $baseDir/bin/rename.sh $aligners_name formatted_designfile.txt
    """

}
process Sort {
    publishDir "${params.outdir}/samtools_sort/", mode: 'link', overwrite: true

    input:
    file( bam_query_file ) from rename_bam_file.collect()

    output:
    file "*_sort*" into sort_bam

    script:
    if (!params.skip_sort){
        println LikeletUtils.print_purple("Samtools sorting the bam files now")
        """    
        bash $baseDir/bin/samtools_sort.sh ${task.cpus}
        """
    } else {
        println LikeletUtils.print_purple("The step of samtools sort is skipped")
        '''
        for bam_file in *.bam
        do
        {
            mv ${bam_file} ${bam_file/.bam/_sort.bam}
        }
        done
        '''
    }
   
}
/*
 * STEP 3-2 - RSeQC analysis
*/
process RSeQC {
    publishDir "${params.outdir}/rseqc" , mode: 'copy', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
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
    /* 
    def strandRule = ''
    if (forward_stranded && !unstranded){
        strandRule = params.singleEnd ? '-d ++,--' : '-d 1++,1--,2+-,2-+'
    } else if (reverse_stranded && !unstranded){
        strandRule = params.singleEnd ? '-d +-,-+' : '-d 1+-,1-+,2++,2--'
    }
    */
    """    
    bash $baseDir/bin/rseqc.sh $bed12 ${task.cpus}
    """
}

process CreateBigWig {
    publishDir "${params.outdir}/rseqc/bigwig", mode: 'link', overwrite: true 

    input:
    file bam from sort_bam.collect()

    output:
    file "*.bigwig" into bigwig_for_genebody

    when:
    !params.skip_rseqc && !params.skip_genebody_coverage 

    script:
    """
    bash $baseDir/bin/create_bigwig.sh ${task.cpus}
    """
}

process GenebodyCoverage {
       publishDir "${params.outdir}/rseqc" , mode: 'link', overwrite: true, 
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
/*
========================================================================================
                            Step 4 Peak Calling
========================================================================================
*/ 
/*
 * STEP 4 - 1  Peak Calling------MetPeak, MACS2, MATK
*/
process Metpeak {
    publishDir "${params.outdir}/peak_calling/metpeak", mode: 'link', overwrite: true

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
    publishDir "${params.outdir}/peak_calling/macs2", mode: 'link', overwrite: true

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
    publishDir "${params.outdir}/peak_calling/MATK", mode: 'link', overwrite: true

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
    matk_jar = baseDir + "/bin/MATK-1.0.jar"
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
    publishDir "${params.outdir}/peak_calling/meyer", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file chromsizesfile from chromsizesfile

    output:
    file "meyer*.bed" into meyer_results
    file "meyer*_normalized.bed" into meyer_nomarlized_bed

    when:
    chromsizesfile != "" && !params.skip_meyer && !params.skip_peakCalling

    shell:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    if( flag_peakCallingbygroup ){
        println LikeletUtils.print_purple("Peak Calling performed by Meyer in group mode")
    }else{
        println LikeletUtils.print_purple("Peak Calling performed by Meyer in independent mode")
    }
    '''
    ln !{baseDir}/bin/meyer.py ./
    bedtools makewindows -g !{chromsizesfile} -w 25 > genome.bin25.bed
    awk '{print $1}' !{chromsizesfile} > chrName.txt
    peak_windows_number=$(wc -l genome.bin25.bed| cut -d " " -f 1)
    mkdir genomebin
    awk '{print "cat genome.bin25.bed | grep "$1" > genomebin/"$1".bin25.bed"}' chrName.txt | xargs -iCMD -P!{task.cpus} bash -c CMD
    bash !{baseDir}/bin/meyer_peakCalling.sh !{formatted_designfile} chrName.txt genomebin/ ${peak_windows_number} !{task.cpus} !{flag_peakCallingbygroup}
    ''' 
}

/*
========================================================================================
                        Step 5 Differential expression analysis
========================================================================================
*/
process Htseq_count{
    publishDir "${params.outdir}/diff_expression/htseq_count", mode: 'link', overwrite: true

    input:
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file gtf

    output:
    file "*input*.count" into htseq_count_input_to_deseq2, htseq_count_input_to_edgeR, htseq_count_input_to_arrange
    file "*ip*.count" into htseq_count_ip_to_arrange

    when:
    !params.skip_expression

    script:
    println LikeletUtils.print_purple("Generate gene expression matrix by htseq-count and Rscript")
    strand_info = params.strand == "no" ? "no" : params.strand == "reverse" ? "reverse" : "yes"
    """
    bash $baseDir/bin/htseq_count.sh $gtf $strand_info ${task.cpus} 
    Rscript $baseDir/bin/get_htseq_matrix.R $formatted_designfile  ${task.cpus} 
    """ 
}

process Deseq2{
    tag "$compare_str"

    publishDir "${params.outdir}/diff_expression/deseq2", mode: 'link', overwrite: true

    input:
    file reads_count_input from htseq_count_input_to_deseq2
    file formatted_designfile from formatted_designfile.collect()
    val compare_str from compareLines_for_Deseq2

    output:
    file "Deseq2*.csv" into deseq2_results
    
    when:
    !params.skip_deseq2 && !params.skip_expression
    
    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by Deseq2 ($compare_str)")
    """
    Rscript $baseDir/bin/DESeq2.R $formatted_designfile $compare_str
    """ 
}

process EdgeR{
    tag "$compare_str"
    publishDir "${params.outdir}/diff_expression/edgeR", mode: 'link', overwrite: true

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
    publishDir "${params.outdir}/diff_expression/cufflinks", mode: 'link', overwrite: true

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
    .into {merged_bed ;mspc_merge_peak_bed; bedtools_merge_peak_bed; bed_for_motif_searching; bed_for_annotation}

process PeakMergeBYMspc {
    publishDir "${params.outdir}/result_arranged/merged_bed", mode: 'link', overwrite: true
    
    input:
    file peak_bed from mspc_merge_peak_bed.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "mspc_group*.bed" into mspc_group_merged_bed
    file "mspc_merged_peaks.bed" into mspc_merged_bed

    when:
    params.mergepeak_mode == "mspc"

    script:
    mspc_directory = baseDir + "/bin/mspc_v3.3"
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    println LikeletUtils.print_purple("Start merge peaks by MSPC")
    """
    bash $baseDir/bin/merge_peaks_by_mspc.sh $mspc_directory $formatted_designfile ${task.cpus} $flag_peakCallingbygroup
    """
}

process PeakMergeBYbedtools {
    publishDir "${params.outdir}/result_arranged/merged_bed", mode: 'link', overwrite: true
    
    input:
    file peak_bed from bedtools_merge_peak_bed.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "bedtools_group*.bed" into bedtools_group_merged_bed
    file "bedtools_merged_peaks.bed" into bedtools_merged_bed

    when:
    params.mergepeak_mode == "bedtools"

    script:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    println LikeletUtils.print_purple("Start merge peaks by bedtools")
    """
    bash $baseDir/bin/merge_peaks_by_bedtools.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup
    """
}
Channel
    .from()
    .concat( mspc_merged_bed, bedtools_merged_bed )
    .into{ merged_bed_all_merged; merged_bed_all_merged_for_diffm6A; merged_bed_part_one }
Channel
    .from()
    .concat( mspc_group_merged_bed, bedtools_group_merged_bed )
    .into{ merged_bed_group_merged_for_diffm6A; merged_bed_group_merged_for_motif; merged_bed_group_merged_for_prediction ; merged_bed_part_two }
    
process PeaksQuantification{
    publishDir "${params.outdir}/result_arranged/quantification", mode: 'link', overwrite: true
    
    input:
    file peak_bed from merged_bed_all_merged.collect()
    file bam_bai_file from sort_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file gtf

    output:
    file "*.{matrix,count}" into quantification_results, quantification_matrix

    when:
    !params.skip_peakCalling

    script:
    matk_jar = baseDir + "/bin/MATK-1.0.jar"
    quantification_mode = params.quantification_mode
    if ( quantification_mode == "bedtools" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by bedtools")
    else if ( quantification_mode == "QNB" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by QNB")
    else if ( quantification_mode == "MATK" )
        println LikeletUtils.print_purple("Generate m6A quantification matrix by MATK")
    """
    # PeaksQuantification by bedtools
    if [ ${quantification_mode} == "bedtools" ]; then 
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${peak_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/bedtools_quantification.R $formatted_designfile bam_stat_summary.txt
    fi

    # PeaksQuantification by QNB
    if [ ${quantification_mode} == "QNB" ]||[ ${quantification_mode} == "MeTDiff" ]; then 
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${peak_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/QNB_quantification.R $formatted_designfile ${task.cpus}
    fi

    # PeaksQuantification by MATK
    if [ ${quantification_mode} == "MATK" ]; then 
        export OMP_NUM_THREADS=1
        bash $baseDir/bin/MATK_quantification.sh $matk_jar $gtf $formatted_designfile ${peak_bed} ${task.cpus}
    fi
    """
}

process diffm6APeak{
    tag "$compare_str"
    publishDir "${params.outdir}/result_arranged/diffm6A", mode: 'link', overwrite: true
    
    input:
    file peak_bed from merged_bed_group_merged_for_diffm6A.collect()
    file merged_bed from merged_bed_all_merged_for_diffm6A.collect()
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
    matk_jar = baseDir + "/bin/MATK-1.0.jar"
    quantification_mode = params.quantification_mode
    if ( quantification_mode == "bedtools" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by bedtools")
    else if ( quantification_mode == "QNB" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by QNB")
    else if ( quantification_mode == "MATK" )
        println LikeletUtils.print_purple("Differential m6A analysis is going on by MATK")
    """
    # PeaksQuantification by bedtools
    if [ ${quantification_mode} == "bedtools" ]; then 
        Rscript $baseDir/bin/bedtools_diffm6A.R $formatted_designfile bedtools_quantification.matrix $compare_str
    fi

    # PeaksQuantification by QNB
    if [ ${quantification_mode} == "QNB" ]; then 
        Rscript $baseDir/bin/QNB_diffm6A.R $formatted_designfile $compare_str 
    fi

    # PeaksQuantification by MeTDiff
    if [ ${quantification_mode} == "MeTdiff" ]; then
        Rscript $baseDir/bin/MeTDiff_diffm6A.R $formatted_designfile $merged_bed $compare_str 
    fi

    # PeaksQuantification by MATK
    if [ ${quantification_mode} == "MATK" ]; then 
        export OMP_NUM_THREADS=${task.cpus}
        bash $baseDir/bin/MATK_diffm6A.sh $matk_jar $formatted_designfile $gtf $compare_str 
    fi

    """ 
}

process MotifSearching {
    publishDir "${params.outdir}/result_arranged/motif", mode: 'link', overwrite: true
    
    input:
    file peak_bed from bed_for_motif_searching.collect()
    file group_bed from merged_bed_group_merged_for_motif.collect()
    file formatted_designfile from formatted_designfile.collect()
    file fasta
    file gtf

    output:
    file "*_dreme" into motif_results

    when:
    !params.skip_dreme

    script:
    println LikeletUtils.print_purple("Motif analysis is going on by DREME and Homer")
    """
    
    bash $baseDir/bin/motif_searching.sh $fasta $gtf ${task.cpus}
    """
}
process SingleNucleotidePrediction{
    publishDir "${params.outdir}/result_arranged/m6Aprediction", mode: 'link', overwrite: true
    
    input:
    file group_bed from merged_bed_group_merged_for_prediction.collect()
    file formatted_designfile from formatted_designfile.collect()
    file fasta
    file gtf

    output:
    file "*_prediction.txt" into prediction_results

    when:
    !params.skip_m6Aprediction

    script:
    matk_jar = baseDir + "/bin/MATK-1.0.jar"
    sramp_dir = baseDir + "/sramp_simple"
    println LikeletUtils.print_purple("SignleNucleotide Prediction analysis is going on by MATK and SRAMP")
    """
    ln -s ${sramp_dir}/* ./
    bash m6Aprediction.sh $matk_jar $fasta $gtf ${task.cpus}
    """
}

Channel
    .from()
    .concat( bed_for_annotation, merged_bed_part_one, merged_bed_part_two )
    .into { annotate_collection; bed_collect_for_arrange_results}

process BedAnnotated{
    publishDir "${params.outdir}/result_arranged/annotation", mode: 'link', overwrite: true
    
    input:
    file annotate_file from annotate_collection.collect()
    file formatted_designfile from formatted_designfile.collect()
    file fasta
    file gtf

    output:
    file "annotatedby{xy,homer}/*" into annotation_results
    
    when:
    !params.skip_annotation

    script:
    annotated_script_dir = baseDir + "/bin"
    //Skip Peak Calling Tools Setting
    """
    # Annotation Peaks
    mkdir -p annotatedbyxy
    mkdir -p annotatedbyhomer
    ln ${annotated_script_dir}/intersec.pl ./
    ln ${annotated_script_dir}/m6A_annotate_forGTF_xingyang2.pl ./
    bash ${baseDir}/bin/annotation.sh ${fasta} ${gtf} ${task.cpus}  
    mv *annotatedbyhomer.bed annotatedbyhomer/
    """
}

Channel
    .from()
    .concat( quantification_results, motif_results, diffm6A_results, 
        htseq_count_ip_to_arrange, htseq_count_input_to_arrange, 
        annotation_results, prediction_results, bed_collect_for_arrange_results
    )
    .set{ results_arrange }

process AggrangeForM6Aviewer {
    publishDir "${params.outdir}/result_arranged/final_results", mode: 'link', overwrite: true
    
    input:
    file results from results_arrange.collect()
    file formatted_designfile from formatted_designfile.collect()
    val compare_info from compareLines_for_arranged_result.collect()
    
    output:
    file "arranged_results.m6APipe" into final_results

    when:
    true

    script:
    mergepeak_mode = params.mergepeak_mode
    quantification_mode = params.quantification_mode
    """
    echo $compare_info | sed 's/^\\[//g' | sed 's/\\]\$//g' | sed s/[[:space:]]//g > compare_info
    Rscript $baseDir/bin/arranged_results.R $formatted_designfile compare_info $mergepeak_mode $quantification_mode
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
   println LikeletUtils.print_yellow("Oops... Pipeline execution stopped with the following message: ")+print_red(workflow.errorMessage)
}