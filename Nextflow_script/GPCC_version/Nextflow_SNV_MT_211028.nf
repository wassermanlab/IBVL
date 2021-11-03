#!/usr/bin/env nextflow

//Nextflow script to process the data for the Silent Genomes Project
// Developped by S. Correard, in October 2021
// Overview of the pipeline :
// For the common part : Fastq --> Bam (using bwa mem)
// SNV calling : SNV calling (Called with DeepVariant) --> Variant merging (using GLNexus)  
//	Extracted MT reads from Bam --> MT reads reverted to fastq --> Reads aligned to the MT ref and the shifted MT ref --> MT Variants calling (using Mutect2)
// Common part again : Variant annotation using VEP

//The version will be defined by the date when the process was started (As it will be a several days process, start day and end day may be different) - format : AAAAMMDD
params.version="20211021"

// Define the path
//Path to fastq folder (input)
params.fastq="/mnt/scratch/SILENT/Act3/GSC_data/Dry_run/*.{R1,R2}.fastq.gz"
fastq_file = file(params.fastq)
//Path to Processed_Individual (output and input)
params.outdir_ind = '/mnt/scratch/SILENT/Act3/Processed/Individual/'
//Path to Processed_Population (output)
params.outdir_pop = '/mnt/scratch/SILENT/Act3/Processed/Population/'


// The reference genome file must be in the scratch space to be loaded by singularity (singularity cannot load a file in the cvmfs)
// Get the genome files
params.ref_genome="/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ref_genome_file = file(params.ref_genome)

//For the bwa mem alignemnt step, the ref genome can be located on cvmfs as it is not a singularity process
// Better to use the one on cvmfs as it is indexed by bwa already
params.ref_genome_cvmfs="/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ref_genome_cvmfs_file = file(params.ref_genome_cvmfs)


//////////////////////////////////////////////////////////////
// Pre-alignment QC
// Fastqc

process fastqc {
	memory '4G'

	publishDir "$params.outdir_ind/${version}/QC/Fastqc/", mode: 'copy'

	input :
	file fastq_file
	val version from params.version
	path outdir_ind from params.outdir_ind

	 output :
	file '*.html' into fastqc_output

	script:
	"""
	# Unload bcchr, and load cvmfs
	# unload_bcchr
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	# load cvmfs
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	
	module load StdEnv/2020
	module load fastqc

	fastqc ${fastq_file} 
	"""
}



//////////////////////////////////////////////////////////////
//Alignment, bam generation

Channel
    .fromFilePairs(params.fastq)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { samples_ch }

// In order to view the files added to the channel, it is possible to add .view() but it HAS TO BE the line after fromFilePairs (Second line)

//
// Step Alignment 1. fastq alignment with bwa mem 
//	Sort with samtools
//

process sorted_bam_files {
	memory '40G'
	cpus 8

	publishDir "$params.outdir_ind/${version}/BAM/", mode: 'copy'

	input :
	ref_genome_cvmfs_file
	set sampleId, file(reads) from samples_ch
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*.bam' optional true into sorted_bam_files
	file '*.bai' optional true into sorted_bam_files_index

	script:
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${ref_genome_cvmfs_file} ${reads} | samtools view -Sb | samtools sort -o ${sampleId}_sorted.bam
	
	samtools index ${sampleId}_sorted.bam
	"""
}

sorted_bam_files.into {
	sorted_bam_files1
	sorted_bam_files2
	sorted_bam_files4
	sorted_bam_files5
	sorted_bam_files6
	sorted_bam_files7
}

sorted_bam_files_index.into {
	sorted_bam_files_index1
	sorted_bam_files_index2
	sorted_bam_files_index4
	sorted_bam_files_index5
	sorted_bam_files_index6
	sorted_bam_files_index7
}

//////////////////////////////////////////////////////////////
//Post-alignment QC

//
// Quality check of the bam file with Mosdepth
//

process Mosdepth {
	memory '4G'

	publishDir "$params.outdir_ind/${version}/QC/Mosdepth/", mode: 'copy'

	input :
	ref_genome_cvmfs_file
	file (bam) from sorted_bam_files7
	file (bai) from sorted_bam_files_index7
	val version from params.version

	output : 
	file '*' optional true into mosdepth

	script :
	"""
	module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/mosdepth-0.3.2.sif \
	mosdepth ${bam.simpleName} ${bam}
	"""
}



//
// Quality check of the bam file with Picard CollectWgsMetrics
// This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses.
//

process Picard_CollectWgsMetrics {
        cpus 8
	publishDir "$params.outdir_ind/${version}/QC/Picard_CollectWgsMetrics/", mode: 'copy'

	input :
	ref_genome_cvmfs_file
	file (bam) from sorted_bam_files6
	file (bai) from sorted_bam_files_index6
	val version from params.version

	output :
	file '*_collect_wgs_metrics.txt' optional true into Picard_CollectWgsMetrics

	script :
	"""
	module load singularity

        singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk	CollectWgsMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_collect_wgs_metrics.txt \
	-R ${ref_genome_cvmfs_file}
	"""
}


//
// Quality check of the bam file with Picard CollectAlignmentSummaryMetrics
// This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters. 
//

process Picard_CollectAlignmentSummaryMetrics {
        cpus 8
	publishDir "$params.outdir_ind/${version}/QC/Picard_CollectAlignmentSummaryMetrics/", mode: 'copy'

	input :
	file (bam) from sorted_bam_files4
	file (bai) from sorted_bam_files_index4
	val version from params.version

	output :
	file '*_Picard_Alignment*' optional true into Picard_CollectAlignmentSummaryMetrics

	script :
	"""
	module load singularity

        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk CollectAlignmentSummaryMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_Picard_Alignment
	"""
}

//
// Quality check of the bam file with Picard QualityScoreDistribution
// This tool is used for determining the overall 'quality' for a library in a given run. To that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases corresponding to those scores. 
// Because R is necessary for the chart, it must be done through Compute Canada
//

process Picard_QualityScoreDistribution {
        cpus 8
	publishDir "$params.outdir_ind/${version}/QC/Picard_QualityScoreDistribution/", mode: 'copy'

	input :
	file (bam) from sorted_bam_files5
	file (bai) from sorted_bam_files_index5
	val version from params.version

	output :
	file '*_qual_score_dist.*' optional true into Picard_QualityScoreDistribution

	script :
	"""
	# Unload bcchr, and load cvmfs
	# unload_bcchr
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	# load cvmfs
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	
	module load StdEnv/2020	
	module load picard
	module load r/4.1.0
	
	java -jar \$EBROOTPICARD/picard.jar \
	QualityScoreDistribution \
	I=${bam} \
	O=${bam.simpleName}_qual_score_dist.txt \
	CHART= ${bam.simpleName}_qual_score_dist.pdf
	"""
}


//
// MultiQC to aggregate all the QC data
//



process multiqc_indiv {
	memory '4G'

	input :
	file '*' from fastqc_output.collect()
	file '*' from Picard_CollectWgsMetrics.collect()
	file '*' from mosdepth.collect()
	file '*' from Picard_CollectAlignmentSummaryMetrics.collect()
	file '*' from Picard_QualityScoreDistribution.collect()
	val version from params.version

	output :
	file '*' optional true into multiqc_indiv

	publishDir "$params.outdir_ind/${version}/QC/multiqc/", mode: 'copy'

	script :
	"""
	module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
	multiqc ${params.outdir_ind}/${version}/QC/
	"""
}






//////////////////////////////////////////////////////////////
//SNV Calling

//
// Step SNV calling 1. Call variants using deepvariant
//


process deepvariant_call {
	memory '40G'
	cpus 8

	input :
	ref_genome_file
	file (bam) from sorted_bam_files1
	file (index) from sorted_bam_files_index1
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	publishDir "$params.outdir_ind/${version}/DeepVariant/", mode: 'copy'

	output :
	file '*.vcf.gz' optional true into deepvariant_call

	script:
	"""
	module load singularity

	echo ${bam.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ -B /mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/ /mnt/common/SILENT/Act3/singularity/deepvariant-1.2.0.sif \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${ref_genome_file} \
	--reads=${bam.simpleName}.bam \
	--regions chr20 \
	--output_gvcf=${bam.simpleName}.g.vcf.gz \
	--output_vcf=${bam.simpleName}.vcf.gz
	"""
}

//
// Step SNV Calling 2. Create a list with all the vcf 
// (!) Need to be done once all calls have been performed
//

process gvcfs_txt {
	input :
	file '*.g.vcf.gz'  from deepvariant_call.collect()
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

 	publishDir "$params.outdir_pop/${version}/", mode: 'copy'

	output :
	file '*.txt' optional true into gvcfs_txt

	script:
	"""
	find $params.outdir_ind/${version}/DeepVariant/ -name "*.g.vcf.gz" > gvcfs.txt
	"""
}

//
// Step SNV Caling 3. GLnexus to combine individual gvcf
//

process GLnexus_cli {
	memory '40G'
	cpus 8

	input :
	file(list_gvcf) from gvcfs_txt
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*.bcf' optional true into GLnexus_cli

	script :
	"""
	module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/glnexus-1.4.1.sif \
	glnexus_cli \
	--config DeepVariant \
	--list ${list_gvcf} > DeepVariant_GLnexus_${version}.bcf
	"""
}

//
// Step SNV Calling 4. Transform the bcf into a vcf
//

process bcf_to_vcf {
	memory '40G'
	cpus 8

	input :
	file(bcf_file) from GLnexus_cli
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*.vcf.gz' into bcf_to_vcf
	file '*.vcf.gz.tbi' into bcf_to_vcf_index

        publishDir "$params.outdir_pop/${version}/", mode: 'copy'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	echo ${bcf_file.simpleName}
	bcftools view ${bcf_file} | bgzip -c > ${bcf_file.simpleName}.vcf.gz

        module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk  --java-options "-Xmx4G" \
       IndexFeatureFile \
        -I ${bcf_file.simpleName}.vcf.gz
	"""
}

	

//////////////////////////////////////////////////////////////////////////////////

//Mitochondrial processes

// Get the genome files for the Mitochondrial data calling (non shifted and shifted)
params.ref_genome_MT="/mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.fasta"
params.ref_genome_MT_shifted="//mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
params.ref_genome_MT_index="//mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.fasta.fai"
params.ref_genome_MT_shifted_index="//mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai"
params.ref_genome_MT_dict="//mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.dict"
params.ref_genome_MT_shifted_dict="//mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict"

ref_genome_MT_file = file(params.ref_genome_MT)
ref_genome_MT_shifted_file = file(params.ref_genome_MT_shifted)
ref_genome_MT_file_index = file(params.ref_genome_MT_index)
ref_genome_MT_shifted_file_index = file(params.ref_genome_MT_shifted_index)
ref_genome_MT_file_dict = file(params.ref_genome_MT_dict)
ref_genome_MT_shifted_file_dict = file(params.ref_genome_MT_shifted_dict)


//
// Step MT 1. index both Ref (shifted by 8000nt and non shifted)
// Fasta downloaded from https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references
//

process MT_Index_Reference {
	input:
	file ref_genome_MT_file
	file ref_genome_MT_shifted_file

	output :
	file '*' into MT_Index_Reference

	script:
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bwa index ${ref_genome_MT_file}
	bwa index ${ref_genome_MT_shifted_file}
	"""
}

MT_Index_Reference.into {
	MT_Index_Reference1
 	MT_Index_Reference2
	MT_Index_Reference3
	MT_Index_Reference4
}

//
// Step MT 2. Extract the MT reads from the bam file
//

process Extract_MT_Read {
	memory '4G'

	input :
	file (bam) from sorted_bam_files2
	file (index) from sorted_bam_files_index2

	output :
	file '*_chrM.bam' into Extract_MT_Read

	script :
	"""
	module load singularity

	echo ${bam.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk PrintReads \
	-L chrM \
	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--read-filter MateUnmappedAndUnmappedReadFilter \
	-I ${bam.simpleName}.bam \
	--read-index ${bam.simpleName}.bam.bai \
	-O ${bam.simpleName}_chrM.bam
	"""
}

process MT_Revert_Sam {
	memory '4G'

	input :
	file(chr_bam) from Extract_MT_Read

	output :
	file '*_chrM_RevertSam.bam' into MT_Revert_Sam

	script :
	"""
	module load singularity
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	echo ${chr_bam.baseName} 

        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk RevertSam \
 	INPUT=${chr_bam.baseName}.bam \
 	OUTPUT_BY_READGROUP=false \
 	OUTPUT=${chr_bam.baseName}_RevertSam.bam \
 	VALIDATION_STRINGENCY=LENIENT \
 	ATTRIBUTE_TO_CLEAR=FT \
 	ATTRIBUTE_TO_CLEAR=CO \
  	SORT_ORDER=queryname \
 	RESTORE_ORIGINAL_QUALITIES=false

	samtools index ${chr_bam.baseName}_RevertSam.bam
	"""
}

process MT_SamtoFastq {
	memory '4G'

	input :
	file(revertSam_bam) from MT_Revert_Sam

	output :
	file '*_chrM_RevertSam.fastq' into MT_SamtoFastq

	script :
	"""
	module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk SamToFastq \
	INPUT=${revertSam_bam.baseName}.bam \
	FASTQ=${revertSam_bam.baseName}.fastq \
	INTERLEAVE=true \
	NON_PF=true
	"""
}

MT_SamtoFastq.into {
  MT_SamtoFastq1
  MT_SamtoFastq2
}

//
//Step MT 3. Align the reads against the ref genome and the shifted ref genome//
//

process align_to_MT {
	memory '4G'

	input :
	path ref_genome_MT_file from MT_Index_Reference1
	file(fastqfromsam) from MT_SamtoFastq1

	output :
	file '*_chrM.bam' into align_to_MT_bam
	file '*_chrM.bam.bai' into align_to_MT_bai

	script:
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	echo ${fastqfromsam.baseName}
	bwa mem -R "@RG\\tID:${fastqfromsam.baseName}\\tSM:${fastqfromsam.baseName}\\tPL:illumina" Homo_sapiens_assembly38.chrM.fasta ${fastqfromsam.baseName}.fastq | samtools view -u -bS | samtools sort > ${fastqfromsam.baseName}_chrM.bam

	samtools index ${fastqfromsam.baseName}_chrM.bam
	"""
}

process align_to_shifted_MT {
	memory '4G'

	input :
	path ref_genome_MT_shifted_file from MT_Index_Reference2
	file(fastqfromsam) from MT_SamtoFastq2

	output :
	file '*_chrM_shifted.bam' into align_to_shifted_MT_bam
	file '*_chrM_shifted.bam.bai' into align_to_shifted_MT_bai

	script:
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bwa mem -R "@RG\\tID:${fastqfromsam.baseName}\\tSM:${fastqfromsam.baseName}\\tPL:illumina" Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta ${fastqfromsam.baseName}.fastq | samtools view -u -bS | samtools sort > ${fastqfromsam.baseName}_chrM_shifted.bam

	samtools index ${fastqfromsam.baseName}_chrM_shifted.bam
	"""
}

//
//Step MT 4. Call the variants against the ref genome and the shifted ref genome
//

process MT_call_variants {
	memory '4G'

	input :
	file ref_genome_MT_file
	file ref_genome_MT_file_index
	file ref_genome_MT_file_dict
	path ref_genome_MT_file from MT_Index_Reference3
	file(bam_MT) from align_to_MT_bam
	file(bai_MT) from align_to_MT_bai

	output :
	file '*_Mutect2.vcf.gz' into MT_call_variants
	file '*_Mutect2.vcf.gz.*' into MT_call_variants_index

	script:
	"""
	module load singularity

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk Mutect2 \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-I ${bam_MT.simpleName}.bam \
	-L chrM \
	--mitochondria-mode \
	--annotation StrandBiasBySample \
	--max-reads-per-alignment-start 75 \
	--max-mnp-distance 0 \
	-O ${bam_MT.simpleName}_Mutect2.vcf.gz
	"""
}

process MT_call_variants_shifted {
	memory '4G'

	input :
	file ref_genome_MT_shifted_file
	file ref_genome_MT_shifted_file_index
	file ref_genome_MT_shifted_file_dict
	path ref_genome_MT_shifted_file from MT_Index_Reference4
	file(shifted_MT_bam) from align_to_shifted_MT_bam
	file(shifted_MT_bai) from align_to_shifted_MT_bai
	
	output :
	file '*_Mutect2.vcf.gz' into MT_call_variants_shifted
	file '*_Mutect2.vcf.gz.*' into MT_call_variants_shifted_index

	script:
	"""
	module load singularity

	echo ${shifted_MT_bam.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk Mutect2 \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-I ${shifted_MT_bam.simpleName}.bam \
	-L chrM \
	--mitochondria-mode \
	--annotation StrandBiasBySample \
	--max-reads-per-alignment-start 75 \
	--max-mnp-distance 0 \
	-O ${shifted_MT_bam.simpleName}_Mutect2.vcf.gz
	"""
}

//
//Step MT 5. Filter Mutect2 calls against the ref genome and the shifted ref genome (Could do the merging before in the future to filter only one file)
//


process MT_Filter_Mutect_Calls {
	memory '4G'

	input :
	file ref_genome_MT_file
	file ref_genome_MT_file_index
	file ref_genome_MT_file_dict
	file vcf_chrM from MT_call_variants
	file vcf_chrM_index from MT_call_variants_index

	output :
	file '*_filtered.vcf.gz' into MT_Filter_Mutect_Calls
	file '*_filtered.vcf.gz.*' into MT_Filter_Mutect_Calls_index

	script :
	"""
	module load singularity

	echo ${vcf_chrM.simpleName}

	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk FilterMutectCalls \
	-V ${vcf_chrM.simpleName}.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.fasta \
	--stats ${vcf_chrM.simpleName}.vcf.gz.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${vcf_chrM.simpleName}_filtered.vcf.gz
	"""
}

process MT_shifted_Filter_Mutect_Calls {
	memory '4G'

	input :
	file ref_genome_MT_shifted_file
	file ref_genome_MT_shifted_file_index
	file ref_genome_MT_shifted_file_dict
	file vcf_chrM_shifted from MT_call_variants_shifted
	file vcf_chrM_shifted_index from MT_call_variants_shifted_index

	output :
	file '*_filtered.vcf.gz' into MT_shifted_Filter_Mutect_Calls
	file '*_filtered.vcf.gz.*' into MT_shifted_Filter_Mutect_Calls_index

	script :
	"""
	module load singularity

	echo ${vcf_chrM_shifted.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk FilterMutectCalls \
	-V ${vcf_chrM_shifted.simpleName}.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	--stats ${vcf_chrM_shifted.simpleName}.vcf.gz.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${vcf_chrM_shifted.simpleName}_filtered.vcf.gz
	"""
}

process MT_LeftAlignAndTrimVariants {
	memory '4G'

	input :
	file ref_genome_MT_file
	file ref_genome_MT_file_index
	file ref_genome_MT_file_dict
	file vcf_fiiltered_chrM from MT_Filter_Mutect_Calls
	file vcf_fiiltered_chrM_index from MT_Filter_Mutect_Calls_index

	output :
	file '*_trimmed.vcf.gz' into MT_LeftAlignAndTrimVariants
	file '*_trimmed.vcf.gz.*' into MT_LeftAlignAndTrimVariants_index

	script :
	"""
	module load singularity

	echo ${vcf_fiiltered_chrM.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${vcf_fiiltered_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_chrM.simpleName}_trimmed.vcf.gz \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
	"""
}

process MT_LeftAlignAndTrimVariants_shifted {
	memory '4G'

	input :
	file ref_genome_MT_shifted_file
	file ref_genome_MT_shifted_file_index
	file ref_genome_MT_shifted_file_dict
	file vcf_fiiltered_chrM_shifted from MT_shifted_Filter_Mutect_Calls
	file vcf_fiiltered_chrM_shifted_index from MT_shifted_Filter_Mutect_Calls_index

	output :
	file '*_trimmed.vcf.gz' into MT_LeftAlignAndTrimVariants_shifted
	file '*_trimmed.vcf.gz.*' into MT_LeftAlignAndTrimVariants_shifted_index
	
	script :
	"""
	module load singularity

	echo ${vcf_fiiltered_chrM_shifted.simpleName}
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-V ${vcf_fiiltered_chrM_shifted.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_chrM_shifted.simpleName}_trimmed.vcf.gz \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
	"""
}

//
//Keep variants in the non-control region
//

process MT_keep_nonCR_variants {
	memory '4G'

	input :
	file ref_genome_MT_file
	file ref_genome_MT_file_index
	file ref_genome_MT_file_dict
	file vcf_fiiltered_trimmed_chrM from MT_LeftAlignAndTrimVariants
	file vcf_fiiltered_trimmed_chrM_index from MT_LeftAlignAndTrimVariants_index

	output :
	file '*_NonControlRegion.vcf.gz' into MT_keep_nonCR_variants

	script :
	"""	
	module load singularity

	echo ${vcf_fiiltered_trimmed_chrM.simpleName}
        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${vcf_fiiltered_trimmed_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_trimmed_chrM.simpleName}_NonControlRegion.vcf.gz \
	-L chrM:576-16024
	"""
}

//
//Keep variants in the control region
//

process MT_keep_CR_variants {
	memory '4G'

	input :
	file ref_genome_MT_shifted_file
	file ref_genome_MT_shifted_file_index
	file ref_genome_MT_shifted_file_dict
	file vcf_fiiltered_trimmed_shifted_chrM from MT_LeftAlignAndTrimVariants_shifted
	file vcf_fiiltered_trimmed_shifted_chrM_index from MT_LeftAlignAndTrimVariants_shifted_index


	output :
	file '*_ControlRegion.vcf.gz' into MT_keep_CR_variants

	script :
	"""
	module load singularity

	echo ${vcf_fiiltered_trimmed_shifted_chrM.simpleName}
        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-V ${vcf_fiiltered_trimmed_shifted_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_trimmed_shifted_chrM.simpleName}_ControlRegion.vcf.gz \
	-L chrM:7455-8576
	"""
}

//
// Shift again the variants to have their real position
// Add the shifted back variants to the previous list of variants
// Sort the variants
// Index the vcf.gz files
//

process MT_shift_variants {
	memory '4G'

	input :
	file vcf_fiiltered_trimmed_CR_region_chrM from MT_keep_CR_variants
	file vcf_fiiltered_trimmed_nonCR_region_chrM from MT_keep_nonCR_variants.collect()
	val version from params.version
	path outdir_ind from params.outdir_ind

	output :
	file '*_merged1_sorted.vcf.gz*' into MT_shift_variants

	publishDir "$params.outdir_ind/${version}/Mutect2/", glob :'*_chrM_Mutect2_filtered_trimmed_sorted.vcf.gz*', mode: 'copy'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
        source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
        conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	echo ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}
	sample_name=\$(echo ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName} | cut -d _ -f 1)
	echo \$sample_name
	gzip -cd ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}.vcf.gz | sed '/^#/d'  | gzip -c > ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader.vcf.gz

	gzip -cd ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}.vcf.gz | grep ^# | gzip -c  > ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_Header.vcf.gz

	gzip -cd ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader.vcf.gz | awk ' \$2>=8001 {print \$1"\t"\$2-8000"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} \$2<=8000 {print \$1"\t"\$2+8569"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} ' | gzip -c  > ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader_ShiftedBack.vcf.gz
	
	cat \${sample_name}_sorted_chrM_RevertSam_chrM_Mutect2_filtered_trimmed_NonControlRegion.vcf.gz ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader_ShiftedBack.vcf.gz > ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1.vcf.gz

	bcftools sort ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1.vcf.gz -O z -o ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1_sorted.vcf.gz
	
	module load singularity

         singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk IndexFeatureFile \
        -I ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1_sorted.vcf.gz
	"""
}


//
//Merge the samples
//

//List the files to merge

process MT_vcfs_txt {
	input :
	file '*_merged1_sorted.vcf.gz' from MT_shift_variants.collect()
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	publishDir "$params.outdir_pop/${version}/", mode: 'copy'

	output :
	file '*.txt' into MT_vcfs_txt

	script:
	"""
	find $params.outdir_ind/${version}/Mutect2/ -name "*_merged1_sorted.vcf.gz" > MT_vcfs.txt
	"""
}

process MT_merge_samples {
	 memory '4G'

	input :
	file (list_MT_vcf) from MT_vcfs_txt
	val version from params.version
	path outdir_pop from params.outdir_pop
	
	output:
	file '*_filtered.vcf.gz' into MT_merge_samples
        file '*_filtered.vcf.gz.*' into MT_merge_samples_index

	publishDir "$params.outdir_pop/${version}/", mode: 'copy'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bcftools merge -l ${list_MT_vcf} -O z -o ${version}_merged_chrM_Mutect2_filtered.vcf.gz

	module load singularity 

         singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk --java-options "-Xmx4G" \
       IndexFeatureFile \
        -I ${version}_merged_chrM_Mutect2_filtered.vcf.gz
	"""
}

////////////////////////////////////////
// For both SNV pop file and Mt pop file

////////////////////////////////////////
// QC on the vcf
// BcfTools stats

process Bcftools_stats {
	memory '4G'
	cpus 8

	input :
	file(vcf_file_MT) from MT_merge_samples
        file(vcf_file_MT_index) from MT_merge_samples_index
	file(vcf_file_SNV) from bcf_to_vcf
        file(vcf_file_SNV_index) from bcf_to_vcf_index
	val version from params.version
	path outdir_pop from params.outdir_pop

	output :
	file '*' into Bcftools_stats

	publishDir "$params.outdir_pop/${version}/QC/Bcftools_stats/", mode: 'copy'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bcftools stats ${vcf_file_MT} > ${vcf_file_MT.simpleName}_bcftools_stat
	bcftools stats ${vcf_file_SNV} > ${vcf_file_SNV.simpleName}_bcftools_stat
	"""
}

//VcfTools TsTv_by_count
// Only performed on SNV as : Error: Polyploidy found, and not supported by vcftools for MT variants
// Do not work on GPCC while it works on Genome Canada, error message : Segmentation fault      (core dumped)

//process Vcftools_TsTv_by_count {
//	memory '8G'
//       cpus 8

//	input :
//	file(vcf_file_SNV) from bcf_to_vcf
//	val version from params.version
//	path outdir_pop from params.outdir_pop

//	output :
//	file '*' into Vcftools_TsTv_by_count

//	publishDir "$params.outdir_pop/${version}/QC/Vcftools_TsTv/", mode: 'copy'

//	script :
//	"""
//	# Unload bcchr, and load cvmfs
//	# unload_bcchr
//	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
//	# load cvmfs
//	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
//	
//	module load StdEnv/2020	
//	module load vcftools

//	echo ${vcf_file_SNV.simpleName}
//	vcftools --gzvcf ${vcf_file_SNV} --TsTv-by-count --out ${vcf_file_SNV.simpleName}_Vcftools_TsTv_count
//	"""
//}

//VcfTools TsTv_by_qual

process Vcftools_TsTv_by_qual {
        memory '4G'

        input :
        file(vcf_file_MT) from MT_merge_samples
        file(vcf_file_SNV) from bcf_to_vcf
        val version from params.version
        path outdir_pop from params.outdir_pop

        output :
        file '*' into Vcftools_TsTv_by_qual

        publishDir "$params.outdir_pop/${version}/QC/Vcftools_TsTv/", mode: 'copy'

        script :
        """
        # Unload bcchr, and load cvmfs
        # unload_bcchr
        source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        # load cvmfs
        source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        module load StdEnv/2020
        module load vcftools

        echo ${vcf_file_MT.simpleName}
	echo ${vcf_file_SNV.simpleName}
        vcftools --gzvcf ${vcf_file_MT} --TsTv-by-qual --out ${vcf_file_MT.simpleName}_Vcftools_TsTv_qual
        vcftools --gzvcf ${vcf_file_SNV} --TsTv-by-qual --out ${vcf_file_SNV.simpleName}_Vcftools_TsTv_qual
	"""
}


//////////////////////////////////////////////////////////////
// Create a file with the list of variants and their annotation
// Same process for SNV and MT variants (using vep, obtaining a tsv)

process annotation_table {
        memory '4G'

        input :
        file(vcf_file_MT) from MT_merge_samples
        file(vcf_file_SNV) from bcf_to_vcf        
        val version from params.version
        path outdir_pop from params.outdir_pop

        output :
        file '*_annotation_tab.tsv' into annotation_table
        file '*_VEP_stats*' into vep_stat

        publishDir "$params.outdir_pop/${version}/", mode: 'copy', pattern : '*tab*'
        publishDir "$params.outdir_pop/${version}/QC/VEP/", mode: 'copy', pattern : '*_VEP_stats*'

        script :
        """
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
        source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
        conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

        echo ${vcf_file_MT}
        echo ${vcf_file_MT.simpleName}

	vep \
        -i ${vcf_file_MT} \
        -o ${vcf_file_MT.simpleName}_${version}_annotation_tab.tsv \
	--offline \
        --cache \
        --dir_cache /mnt/common/DATABASES/REFERENCES/GRCh38/VEP/ \
        --everything \
        --tab \
        --stats_file ${vcf_file_MT.simpleName}_VEP_stats

	vep \
        -i ${vcf_file_SNV} \
        -o ${vcf_file_SNV.simpleName}_${version}_annotation_tab.tsv \
        --offline \
	--cache \
        --dir_cache /mnt/common/DATABASES/REFERENCES/GRCh38/VEP/ \
        --everything \
        --tab \
        --stats_file ${vcf_file_SNV.simpleName}_VEP_stats

        """
}


//////////////////////////////////////////////////////////////
// Create a file with the list of SNV and their frequency in the dataset
// Only works for the SNV frequency as other values are necessary for the MT variants


process SNV_frequency_table {
        memory '4G'

        input :
        file(SNV_vcf) from bcf_to_vcf
	file(SNV_vcf_index) from bcf_to_vcf_index
        val version from params.version
        path outdir_pop from params.outdir_pop

        output :
        file '*' into SNV_frequency_table

        publishDir "$params.outdir_pop/${version}/", mode: 'copy'

        script :
        """
	module load singularity

         echo ${SNV_vcf}
         echo ${SNV_vcf.simpleName}

         singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk --java-options "-Xmx4G" \
       VariantsToTable \
        -V ${SNV_vcf} \
        -O ${SNV_vcf.simpleName}_frequency_table \
        -F CHROM \
        -F POS \
        -F TYPE\
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F AF \
        -F HET \
        -F HOM-REF \
        -F HOM-VAR \
        -F NO-CALL \
        -F MULTI-ALLELIC \
        -F Consequence
        """
}


//////////////////////////////////////////////////////////////
// Create a file with the list of MT and their frequency in the dataset
// Only works for the MT frequency as it is necessary to calculate the VAF (Variant allele fraction or heteroplasmy levels)


process MT_frequency_table {
        memory '4G'

        input :
        file(MT_vcf) from MT_merge_samples
        file(MT_vcf_index) from MT_merge_samples_index
	val version from params.version
        path outdir_pop from params.outdir_pop

        output :
        file '*' into MT_frequency_table

        publishDir "$params.outdir_pop/${version}/", mode: 'copy'

        script :
        """
         echo ${MT_vcf}
         echo ${MT_vcf.simpleName}

	module load singularity

        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk --java-options "-Xmx4G" \
       VariantsToTable \
        -V ${MT_vcf} \
        -O ${MT_vcf.simpleName}_frequency_table \
        -F CHROM \
        -F POS \
        -F TYPE\
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F AF \
        -F HET \
        -F HOM-REF \
        -F HOM-VAR \
        -F NO-CALL \
        -F MULTI-ALLELIC \
        -F Consequence \
        -GF AF
        """
}



//
// MultiQC to aggregate all the QC data
//



process multiqc_pop {
        memory '4G'

        input :
        file '*' from vep_stat.collect()
        file '*' from Vcftools_TsTv_by_qual.collect()
        file '*' from Bcftools_stats.collect()
        val version from params.version
        path outdir_pop from params.outdir_pop

        output :
        file '*' into multiqc_pop

        publishDir "$params.outdir_pop/${version}/QC/multiQC/", mode: 'copy'

        script :
        """
        module load singularity

        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
        multiqc $params.outdir_pop/${version}/QC/
        """
}
