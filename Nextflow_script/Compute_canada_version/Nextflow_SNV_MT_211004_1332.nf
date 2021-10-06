#!/usr/bin/env nextflow

//Nextflow script to process the data for the Silent Genomes Project
// Developped by S. Correard, in October 2021
// Overview of the pipeline :
// For the common part : Fastq --> Bam (using bwa mem)
// SNV calling : SNV calling (Called with DeepVariant) --> Variant merging (using GLNexus)  
//	Extracted MT reads from Bam --> MT reads reverted to fastq --> Reads aligned to the MT ref and the shifted MT ref --> MT Variants calling (using Mutect2)
// Common part again : Variant annotation using VEP

//The version will be defined by the date when the process was started (As it will be a several days process, start day and end day may be different) - format : AAAAMMDD
params.version="20211005"

// Define the path
//Path to fastq folder (input)
params.fastq="/home/correard/scratch/big_files/GSC/Batch_1/FASTQ/*_{1,2}.fastq.gz"
//Path to Processed_Individual (output and input)
params.outdir_ind = '/home/correard/scratch/SnakeMake_VariantCalling/Processed/Processed_Individuals/'
//Path to Processed_Population (output)
params.outdir_pop = '/home/correard/scratch/SnakeMake_VariantCalling/Processed/Processed_Population/'


// The reference genome file must be in the scratch space to be loaded by singularity (singularity cannot load a file in the cvmfs)
// Get the genome files
params.ref_genome="/home/correard/scratch/SnakeMake_VariantCalling/snakemake_SNV_full/Homo_sapiens.GRCh37.fa"
ref_genome_file = file(params.ref_genome)


//////////////////////////////////////////////////////////////
// Data organization - Create the appropriate folders

process folder_creation {
	input :
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	script :
	"""
	mkdir -p ${outdir_ind}/BAM/${version}/
	mkdir -p ${outdir_ind}/DeepVariant/${version}/
	mkdir -p ${outdir_ind}/Mutect2/${version}/
	mkdir -p ${outdir_pop}/${version}/
	"""
}


//////////////////////////////////////////////////////////////
//Alignment, bam generation

Channel
    .fromFilePairs(params.fastq)
    .view()
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { samples_ch }

//
// Step Alignment 1. fastq alignment with bwa mem 
//	Sort with samtools
//

process sorted_bam_files {
	memory '40G'
	cpus 24
	time '224h'

	publishDir '${outdir_ind}/BAM/${version}/', mode: 'copy'

	input :
	ref_genome_file
	set sampleId, file(reads) from samples_ch
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*' into sorted_bam_files

	script:
	"""
	bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${ref_genome_file} ${reads} | samtools view -Sb - | samtools sort -o ${sampleId}_sorted.bam
	
	samtools index ${sampleId}_sorted.bam
	"""
}


Channel
    .fromFilePairs('${outdir_ind}/BAM/${version}/*.{bam,bai}') { file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .set { samples_bam_ch }
//	.view()

samples_bam_ch.into {
	samples_bam_ch1
	samples_bam_ch2
	samples_bam_ch3
}


//////////////////////////////////////////////////////////////
//SNV Calling

//
// Step SNV calling 1. Call variants using deepvariant
//


process deepvariant_call {
	memory '40G'
	cpus 24
	time '224h'

	input :
	ref_genome_file
	set sample_bam_id, file(ordered_bam) from samples_bam_ch1
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	publishDir '${outdir_ind}/DeepVariant/${version}/', mode: 'copy'

	output :
	file '*.vcf.gz' into deepvariant_call

	script:
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/snakemake-deepvariant/deepvariant.sif \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${ref_genome_file} \
	--reads=${sample_bam_id}.bam \
	--regions "20" \
	--output_gvcf=${sample_bam_id}.g.vcf.gz \
	--output_vcf=${sample_bam_id}.vcf.gz
	"""
}

//
// Step SNV Calling 2. Create a list with all the vcf 
// (!) Need to be done once all calls have been performed
//

process gvcfs_txt {
	input :
	file '*.g.vcf.gz'  from deepvariant_call.toList()
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	publishDir 	publishDir '${outdir_pop}/${version}/', mode: 'copy'

	output :
	file '*.txt' into gvcfs_txt

	script:
	"""
	find ${outdir_ind}/DeepVariant/ -name "*.g.vcf.gz" > gvcfs.txt
	"""
}

//
// Step SNV Caling 3. GLnexus to combine individual gvcf
//

process GLnexus_cli {
	memory '40G'
	cpus 24

	input :
	file(list_gvcf) from gvcfs_txt
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*.bcf' into GLnexus_cli

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/snakemake-GLnexus/GLnexus.sif \
	glnexus_cli \
	--config DeepVariant \
	--list ${list_gvcf} > ${version}.bcf
	"""
}

//
// Step SNV Calling 4. Transform the bcf into a vcf
//

process bcf_to_vcf {
	memory '40G'
	cpus 24

	input :
	file(bcf_file) from GLnexus_cli
	val version from params.version
	path outdir_ind from params.outdir_ind
	path outdir_pop from params.outdir_pop

	output :
	file '*.vcf.gz' into bcf_to_vcf

	script :
	"""
	bcftools view ${bcf_file} | bgzip -c > ${version}.vcf.gz
	"""
}

	

//////////////////////////////////////////////////////////////////////////////////

//Mitochondrial processes

// Get the genome files for the Mitochondrial data calling (non shifted and shifted)
params.ref_genome="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.fasta"
params.ref_genome_shifted="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
params.ref_genome_index="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.fasta.fai"
params.ref_genome_shifted_index="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai"
params.ref_genome_dict="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.dict"
params.ref_genome_shifted_dict="/home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict"

ref_genome_file = file(params.ref_genome)
ref_genome_shifted_file = file(params.ref_genome_shifted)
ref_genome_file_index = file(params.ref_genome_index)
ref_genome_shifted_file_index = file(params.ref_genome_shifted_index)
ref_genome_file_dict = file(params.ref_genome_dict)
ref_genome_shifted_file_dict = file(params.ref_genome_shifted_dict)


//
// Step MT 1. index both Ref (shifted by 8000nt and non shifted)
// Fasta downloaded from https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references
//

process MT_Index_Reference {
	input:
	file ref_genome_file
	file ref_genome_shifted_file

	output :
	file '*' into MT_Index_Reference

	script:
	"""
	bwa index ${ref_genome_file}
	bwa index ${ref_genome_shifted_file}
	"""
}

MT_Index_Reference.into {
	MT_Index_Reference1
 	MT_Index_Reference2
	MT_Index_Reference3
	MT_Index_Reference4
	MT_Index_Reference5
	MT_Index_Reference6
	MT_Index_Reference7
	MT_Index_Reference8
}

//
// Step MT 2. Extract the MT reads from the bam file
//

process Extract_MT_Read {
	memory '4G'

	input :
	set sample_bam_id, file(ordered_bam) from samples_bam_ch2

	output :
	file '*_chrM.bam' into Extract_MT_Read

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk --java-options "-Xmx4G" \
	PrintReads \
	-L MT \
	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--read-filter MateUnmappedAndUnmappedReadFilter \
	-I ${ordered_bam} \
	--read-index ${sample_id}.bam.bai \
	-O ${sample_id}_chrM.bam
	"""
}

process MT_Revert_Sam {
	memory '4G'

	input :
	path bam_chrM from Extract_MT_Read
	set sample_bam_id, file(ordered_bam) from samples_bam_ch3

	output :
	file '*_chrM_RevertSam.bam' into MT_Revert_Sam

	script :
	"""
	echo ${sample_id} 
	java -jar $EBROOTPICARD/picard.jar \
	RevertSam \
 	INPUT=${sample_id}_chrM.bam \
 	OUTPUT_BY_READGROUP=false \
 	OUTPUT=${sample_id}_chrM_RevertSam.bam \
 	VALIDATION_STRINGENCY=LENIENT \
 	ATTRIBUTE_TO_CLEAR=FT \
 	ATTRIBUTE_TO_CLEAR=CO \
  	SORT_ORDER=queryname \
 	RESTORE_ORIGINAL_QUALITIES=false

	samtools index ${sample_id}_chrM_RevertSam.bam
	"""
}

MT_Revert_Sam.into {
  MT_Revert_Sam1
  MT_Revert_Sam2
}


process MT_SamtoFastq {
	memory '4G'

	input :
	path bam_reverted from MT_Revert_Sam2

	output :
	file '*_chrM_RevertSam.fastq' into MT_SamtoFastq

	script :
	"""
	java -jar $EBROOTPICARD/picard.jar \
	SamtoFastq \
	INPUT=${sample_id}_chrM_RevertSam.bam \
	FASTQ=${sample_id}_chrM_RevertSam.fastq \
	INTERLEAVE=true \
	NON_PF=true
	"""
}

MT_SamtoFastq.into {
  MT_SamtoFastq1
  MT_SamtoFastq2
}

//
//Stap MT 3. Align the reads against the ref genome and the shifted ref genome//
//

process align_to_MT {
	memory '4G'

	input :
	path fastq from MT_SaMToFastq1
	path ref_genome_file MT_Index_Reference1

	output :
	file '*' into align_to_MT

	script:
	 """
	bwa mem -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" Homo_sapiens_assembly38.chrM.fasta ${sample_id}_chrM_RevertSam.fastq | samtools view -u -bS | samtools sort > ${sample_id}_chrM_RevertSam_chrM.bam

	samtools index ${sample_id}_chrM_RevertSam_chrM.bam
	"""
}

process align_to_shifted_MT {
	memory '4G'

	input :
	path fastq from MT_SaMToFastq2
	path ref_genome_shifted_file from MT_Index_Reference2

	output :
	file '*' into align_to_shifted_MT

	script:
	"""
	bwa mem -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta ${sample_id}_chrM_RevertSam.fastq | samtools view -u -bS | samtools sort > ${sample_id}_chrM_RevertSam_chrM_shifted.bam

	samtools index ${sample_id}_chrM_RevertSam_chrM_shifted.bam
	"""
}

//
//Step MT 4. Call the variants against the ref genome and the shifted ref genome
//

process MT_call_variants {
	memory '4G'

	input :
	file ref_genome_file
	file ref_genome_file_index
	file ref_genome_file_dict
	path ref_genome_file from MT_Index_Reference2
	path bam_chrM from align_to_MT

	output :
	file '*' into MT_call_variants

	script:
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk Mutect2 \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-I ${sample_id}_chrM_RevertSam_chrM.bam \
	-L chrM \
	--mitochondria-mode \
	--annotation StrandBiasBySample \
	--max-reads-per-alignment-start 75 \
	--max-mnp-distance 0 \
	-O ${sample_id}_chrM_RevertSam_chrM_Mutect2.vcf.gz
	"""
}

process MT_call_variants_shifted {
	memory '4G'

	input :
	file ref_genome_shifted_file
	file ref_genome_shifted_file_index
	file ref_genome_shifted_file_dict
	path ref_genome_shifted_file from MT_Index_Reference4
	path bam_chrM_shifted from align_to_shifted_MT

	output :
	file '*' into MT_call_variants_shifted

	script:
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk Mutect2 \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-I ${sample_id}_chrM_RevertSam_chrM_shifted.bam \
	-L chrM \
	--mitochondria-mode \
	--annotation StrandBiasBySample \
	--max-reads-per-alignment-start 75 \
	--max-mnp-distance 0 \
	-O ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2.vcf.gz
	"""
}

//
//Step MT 5. Filter Mutect2 calls against the ref genome and the shifted ref genome (Could do the merging before in the future to filter only one file)
//


process MT_Filter_Mutect_Calls {
	memory '4G'

	input :
	file ref_genome_file
	file ref_genome_file_index
	file ref_genome_file_dict
	file vcf_chrM from MT_call_variants
	
	output :
	file '*' into MT_Filter_Mutect_Calls

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk FilterMutectCalls \
	-V ${sample_id}_chrM_RevertSam_chrM_Mutect2.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.fasta \
	--stats ${sample_id}_chrM_RevertSam_chrM_Mutect2.vcf.gz.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered.vcf.gz
	"""
}

process MT_shifted_Filter_Mutect_Calls {
	memory '4G'

	input :
	file ref_genome_shifted_file
	file ref_genome_shifted_file_index
	file ref_genome_shifted_file_dict
	file vcf_chrM from MT_call_variants_shifted

	output :
	file '*' into MT_shifted_Filter_Mutect_Calls

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk FilterMutectCalls \
	-V ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	--stats ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2.vcf.gz.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered.vcf.gz
	"""
}

process MT_LeftAlignAndTrimVariants {
	memory '4G'

	input :
	file ref_genome_file
	file ref_genome_file_index
	file ref_genome_file_dict
	file vcf_fiiltered_chrM from MT_Filter_Mutect_Calls

	output :
	file '*' into MT_LeftAlignAndTrimVariants
	
	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered.vcf.gz \
	-O ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered_trimmed.vcf.gz \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
	"""
}

process MT_LeftAlignAndTrimVariants_shifted {
	memory '4G'

	input :
	file ref_genome_shifted_file
	file ref_genome_shifted_file_index
	file ref_genome_shifted_file_dict
	file vcf_fiiltered_chrM_shifted from MT_shifted_Filter_Mutect_Calls

	output :
	file '*' into MT_LeftAlignAndTrimVariants_shifted

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-V ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered.vcf.gz \
	-O ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed.vcf.gz \
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
	file ref_genome_file
	file ref_genome_file_index
	file ref_genome_file_dict
	file vcf_fiiltered_trimmed_chrM from MT_LeftAlignAndTrimVariants

	output :
	file '*' into MT_keep_nonCR_variants

	script :
	"""	
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered_trimmed.vcf.gz \
	-O ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered_trimmed_NonControlRegion.vcf.gz \
	-L chrM:576-16024
	"""
}

//
//Keep variants in the control region
//

process MT_keep_CR_variants {
	memory '4G'

	input :
	file ref_genome_shifted_file
	file ref_genome_shifted_file_index
	file ref_genome_shifted_file_dict
	file vcf_fiiltered_trimmed_shifted_chrM from MT_LeftAlignAndTrimVariants_shifted

	output :
	file '*' into MT_keep_CR_variants

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-V ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed.vcf.gz \
	-O ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion.vcf.gz \
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
	file vcf_fiiltered_trimmed_nonCR_region_chrM from MT_keep_nonCR_variants

	output :
	file '*' into MT_shift_variants

	publishDir '${outdir_ind}/Mutect2/${version}/', glob :'*_chrM_Mutect2_filtered_trimmed_sorted.vcf.gz', mode: 'copy'

	script :
	"""
	gzip -cd ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion.vcf.gz | sed '/^#/d'  | gzip -c > ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion_NoHeader.vcf.gz

	gzip -cd ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion.vcf.gz | grep ^# | gzip -c  > ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion_Header.vcf.gz

	gzip -cd ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion_NoHeader.vcf.gz | awk ' \$2>=8001 {print \$1"\t"\$2-8000"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} \$2<=8000 {print \$1"\t"\$2+8569"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} ' | gzip -c  > ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion_NoHeader_ShiftedBack.vcf.gz
	
	cat ${sample_id}_chrM_RevertSam_chrM_Mutect2_filtered_trimmed_NonControlRegion.vcf.gz ${sample_id}_chrM_RevertSam_chrM_shifted_Mutect2_filtered_trimmed_ControlRegion_NoHeader_ShiftedBack.vcf.gz > ${sample_id}_chrM_Mutect2_filtered_trimmed.vcf.gz

	bcftools sort ${sample_id}_chrM_Mutect2_filtered_trimmed.vcf.gz -O z -o ${sample_id}_chrM_Mutect2_filtered_trimmed_sorted.vcf.gz

	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
   gatk IndexFeatureFile \
   -I ${sample_id}_chrM_Mutect2_filtered_trimmed_sorted.vcf.gz
	"""
}

//
//Merge the samples
//

process MT_merge_samples {
	 memory '4G'

	input :
	file vcf_individuals from MT_shift_variants
	val version from params.version
	path outdir_pop from params.outdir_pop

	output:
	file '*' into MT_merge_samples

	publishDir '${outdir_pop}/${version}/', mode: 'copy'

	script :
	"""
	bcftools merge *_chrM_Mutect2_filtered_trimmed_sorted.vcf.gz -O z -o merged_chrM_Mutect2_filtered.vcf.gz
	"""
}

////////////////////////////////////////
// For both SNV pop file and Mt pop file
// Annotate the vcf with VEP
// Conpress the annotated vcf
//


process annotate_vcf {
	memory '4G'

	input :
	file(vcf_file) from MT_merge_samples
	file(vcf_file) from bcf_to_vcf
	val version from params.version
	path outdir_pop from params.outdir_pop

	output :
	file '*' into annotate_vcf

	publishDir '${outdir_pop}/${version}/', mode: 'copy'

	script :
	"""
	singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/snakemake-VEP/vep.sif \
	vep \
	-i ${vcf_file} \
	-o ${vcf_file}_${version}_annotated.vcf \
	--cache \
	--dir_cache /home/correard/scratch/SnakeMake_VariantCalling/snakemake-VEP/vep_cache \
	--everything \
	--vcf

	bgzip -c ${vcf_file}_${version}_annotated.vcf > ${vcf_file}_${version}_annotated.vcf.gz
	"""
}


