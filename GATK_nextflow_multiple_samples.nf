//GATK Nextflow pipeline for multiple samples
//Dr Chris J Smith 
//09/01/2025

/*
 * Pipeline parameters
 */

// Execution environment setup
params.scratchDir = "/path/to/working/directory"
params.gatkDir = "/path/to/PublicDataSets/GATKbundle/hg38/v0"
params.outdir_bam = "/path/to/working/directory/output/bam"
params.outdir_vcf = "/path/to/working/directory/output/vcf"

// Primary input for FASTQ files
params.inputDir = "${params.scratchDir}/input"

// Accessory files
params.bed = "${params.scratchDir}/twist_bed_plus_200.bed"
params.ref = "${params.gatkDir}/Homo_sapiens_assembly38.fasta"
params.snp = "${params.gatkDir}/Homo_sapiens_assembly38.dbsnp138.vcf"
params.indel = "${params.gatkDir}/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.indel2 = "${params.gatkDir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

// Create channel from input fastq files
Channel
    .fromFilePairs("${params.inputDir}/*_{1,2}.fastq.gz")
    .map { filename, files -> 
        def sampleName = filename[0..5]  // Extract first 6 chars/change as needed
        tuple(
            files[0],                    // fastq1
            files[1],                    // fastq2
            sampleName,                  
            "READGROUPID",              //Enter RG information here
            "ILLUMINA",                  
            "PLATFORM_UNIT",       
            "LIBRARY"          
        )
    }
    .set { fastq_ch }

// Processes

/*
 * Generate BAM file from FASTQ files
 */

process FastqToSam {

    //Loads GATK before running the script
    beforeScript "module load gatk"
    
    //Input tuple for read group information
    input:
       tuple path(fastq1), 
             path(fastq2), 
             val(sampleName), 
             val(readGroupID), 
             val(platform), 
             val(platformUnit), 
             val(library)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.unaligned.bam"), emit: bam
    
    // Set the command to be executed
    script:
    """
    gatk FastqToSam \
        -F1 ${fastq1} \
        -F2 ${fastq2} \
        -O ${sampleName}.unaligned.bam \
        -RG ${readGroupID} \
        -LB ${library} \
        -PL ${platform} \
        -PU ${platformUnit} \
        -SM ${sampleName}
    """
}

/*
 * Mark Illumina adaptors
 */
process MarkIlluminaAdapters {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName), 
              path(unaligned_bam)
    
    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.added.markilluminaadapters.bam"), emit: bam

    script:
    """
    gatk MarkIlluminaAdapters \
        -I ${unaligned_bam} \
        -O ${sampleName}.added.markilluminaadapters.bam \
        -M ${sampleName}.markilluminaadapters.metrics.txt
    """
}

/*
 * SamtoFastq
 */
process SamtoFastq {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName), 
              path(adapters_bam)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.fq")

      
    script:
    """
    gatk SamToFastq \
        -I ${adapters_bam} \
        -F ${sampleName}.fq \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2 \
        --INTERLEAVE true \
        --NON_PF true;
    """
}

/*
 * Align to reference genome
 */
process bwa_mem {
    //Loads BWA before running the script
    beforeScript "module load bwa"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName), 
              path(input_fq)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.sam")

   script:
    """
    bwa mem -M -t 2 \
    -p ${params.ref} ${input_fq} > ${sampleName}.sam
    """
}

// /*
//  * Merge BAM files
//  */
process MergeBamAlignment {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(unaligned_bam), 
              path(aligned_bam)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.mergebamalignment.bam"), emit: bam

       
    script:

    """
    gatk MergeBamAlignment \
    -R ${params.ref} \
    --UNMAPPED_BAM ${sampleName}.unaligned.bam \
    --ALIGNED_BAM ${sampleName}.sam \
    -O ${sampleName}.mergebamalignment.bam \
    --CREATE_INDEX true \
    --CLIP_ADAPTERS false \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ATTRIBUTES_TO_RETAIN XS;
    """
}

// /*
//  * Mark Duplicates
//  */
process MarkDuplicates {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(mergebamalignment_bam)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.markduplicates.bam"), emit: bam
     
    script:
    """
    gatk MarkDuplicates \
    -I ${sampleName}.mergebamalignment.bam \
    -O ${sampleName}.markduplicates.bam \
    -M ${sampleName}.markduplicates.metrics.txt;
    """
}

// /*
//  * BaseRecalibrator
//  */
process BaseRecalibrator {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(markduplicates_bam)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.markduplicates.table")
     
    script:
    """
    gatk BaseRecalibrator \
    -R ${params.ref} \
    -I ${sampleName}.markduplicates.bam \
    --known-sites ${params.snp} \
    --known-sites ${params.indel} \
    --known-sites ${params.indel} \
    -O ${sampleName}.markduplicates.table;
    """
}

// /*
//  * ApplyBQSR
//  */
process ApplyBQSR {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(markduplicates_bam),
              path(markduplicates_table)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.recal.bam"), 
              path("${sampleName}.recal.bai")
     
     // Copy output files to the output bam folder
    publishDir params.outdir_bam, mode: 'copy'

    script:
    """
    gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${sampleName}.markduplicates.bam \
    --bqsr-recal-file ${sampleName}.markduplicates.table \
    -O ${sampleName}.recal.bam;
    """
}

// /*
//  * HaplotypeCaller
//  */
process HaplotypeCaller {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(recal_bam),
              path(recal_bai)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.raw.vcf")
     
    script:
    """
    gatk HaplotypeCaller \
    -R ${params.ref} \
    -I ${sampleName}.recal.bam \
    -L ${params.bed} \
    --dbsnp ${params.snp} \
    -O ${sampleName}.raw.vcf;
    """
}

// /*
//  * SNP selection
//  */
process SNPSelection {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(raw_vcf)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.snps.vcf")
     
    script:
    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${sampleName}.raw.vcf \
    -select-type SNP \
    -O ${sampleName}.snps.vcf;
    """
}

// /*
//  * SNP filtering
//  */
process SNPFiltering {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(snps_vcf)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.filtered.snps.vcf")
     
    script:
    """
    gatk VariantFiltration \
    -R ${params.ref} \
    -V ${sampleName}.snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter"  \
    -O ${sampleName}.filtered.snps.vcf;
    """
}

// /*
//  * INDEL selection
//  */
process IndelSelection {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(raw_vcf)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.indels.vcf")
     
    script:
    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${sampleName}.raw.vcf \
    -select-type INDEL \
    -O ${sampleName}.indels.vcf;
    """
}

// /*
//  * Indel filtering
//  */
process IndelFiltering {
    beforeScript "module load gatk"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(indel_vcf)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.filtered.indels.vcf")
     
    script:
    """
    gatk VariantFiltration \
    -R ${params.ref} \
    -V ${sampleName}.indels.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"\
    --filter-name "my_indel_filter"  \
    -O ${sampleName}.filtered.indels.vcf;
    """
}

// /*
//  * VCF prosessing
//  */
process VCFProcessing {
        
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(snp_filtered_vcf),
              path(indel_filtered_vcf)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}.filtered.merged.vcf.gz")
     
    // Copy output files to the output vcf folder
    publishDir params.outdir_vcf, mode: 'copy'
    
    script:
    """
    #!/bin/bash
    
    module load bcftools
    module load htslib

    bgzip ${sampleName}.filtered.snps.vcf
    bgzip ${sampleName}.filtered.indels.vcf
    tabix ${sampleName}.filtered.snps.vcf.gz
    tabix ${sampleName}.filtered.indels.vcf.gz

    bcftools concat -a ${sampleName}.filtered.snps.vcf.gz \
                       ${sampleName}.filtered.indels.vcf.gz \
                       -o ${sampleName}.filtered.merged.vcf

    bgzip ${sampleName}.filtered.merged.vcf
    tabix ${sampleName}.filtered.merged.vcf.gz
    """
}

workflow {

    // Map the tuple to the individual parameters for FastqToSam
    FastqToSam(fastq_ch)

    //Run GATK MarkIlluminaAdapters
    MarkIlluminaAdapters(FastqToSam.out)

   //Run GATK SamtoFastq
    SamtoFastq(MarkIlluminaAdapters.out)

    //Run BWA mem
    bwa_mem(SamtoFastq.out)

    //Run GATK MergeBamAlignment
    MergeBamAlignment(FastqToSam.out.join(bwa_mem.out)) //Join the outputs to input unaligned and aligned bams

    //Run GATK MarkDuplicates
    MarkDuplicates(MergeBamAlignment.out)

    //Run GATK BaseRecalibrator
    BaseRecalibrator(MarkDuplicates.out)

    //Run GATK ApplyBQSR
    ApplyBQSR(MarkDuplicates.out.join(BaseRecalibrator.out)) //Join the outputs to input bam and table

    //Run GATK HaplotypeCaller
    HaplotypeCaller(ApplyBQSR.out)

    //Run GATK SNPSelection
    SNPSelection(HaplotypeCaller.out)

    //Run GATK SNPFiltering
    SNPFiltering(SNPSelection.out)

    //Run GATK IndelSelection
    IndelSelection(HaplotypeCaller.out)

    //Run GATK IndelFiltering
    IndelFiltering(IndelSelection.out)

    //Run GATK VCFProcessing
    VCFProcessing(SNPFiltering.out.join(IndelFiltering.out)) //Join the outputs to input snp and indel vcf files

}


