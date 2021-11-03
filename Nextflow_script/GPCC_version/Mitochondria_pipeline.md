
# Mitochondria_pipeline file

This readme is for everyone to comment on things to add / modify in the script

This Readme only concerns the Mitochondrial part of the pipeline. Some parts are common with the SNV pipeline and therefore are described in the [alignment_SNV_pipeline](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/Alignment_SNV_pipeline.md) document.

## References downloading and indexing

It is necessary to download the reference genome as well as the shifted reference genomes

Fasta downloaded from https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references

```
	bwa index ${ref_genome_MT_file}
	bwa index ${ref_genome_MT_shifted_file}
 ```

## MT variant calling
  
  Decided to use GATK and follow the broad guidelines : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  It requires many step as reads mapping to the MT genome are extracted, and then re-aligned to a MT genome and a shifted MT genome to adress the circularity of the genome (And the reads mapping on the artificial breakpoint in the linear reference genome).
  
  The steps are described here : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  The general idea was kept while some step are slightly different


 
 ### Extract MT reads
 
 May be necessary to also extract the reads that maps to NuMTs
 
 ```
 	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk PrintReads \
	-L chrM \
	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--read-filter MateUnmappedAndUnmappedReadFilter \
	-I ${bam.simpleName}.bam \
	--read-index ${bam.simpleName}.bam.bai \
	-O ${bam.simpleName}_chrM.bam
  ```
  
  ### Remove all the position information from the bam while keeping other information
  
  ```
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
```

### Change bam file to fastq

```
	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk SamToFastq \
	INPUT=${revertSam_bam.baseName}.bam \
	FASTQ=${revertSam_bam.baseName}.fastq \
	INTERLEAVE=true \
	NON_PF=true
 ```
 
 ### Align Fastq to MT reference genome (also run using the shifted MT reference genome)
 
 ```
 	bwa mem -R "@RG\\tID:${fastqfromsam.baseName}\\tSM:${fastqfromsam.baseName}\\tPL:illumina" Homo_sapiens_assembly38.chrM.fasta ${fastqfromsam.baseName}.fastq | samtools view -u -bS | samtools sort > ${fastqfromsam.baseName}_chrM.bam
  
	samtools index ${fastqfromsam.baseName}_chrM.bam
 ```
 
 ### Call MT variants using Mutect2 (also run using the shifted MT reference genome)
 
 ```
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
  ```
  
  ### Filter the calls (also run using the shifted MT reference genome)
  
 ```
 	singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk FilterMutectCalls \
	-V ${vcf_chrM.simpleName}.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.fasta \
	--stats ${vcf_chrM.simpleName}.vcf.gz.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${vcf_chrM.simpleName}_filtered.vcf.gz
  ```
  
  ### LeftAlignAndTrimVariants (also run using the shifted MT reference genome)
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${vcf_fiiltered_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_chrM.simpleName}_trimmed.vcf.gz \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
 ```
 
 ### Keep variants located in the non CR region
 
 ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${vcf_fiiltered_trimmed_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_trimmed_chrM.simpleName}_NonControlRegion.vcf.gz \
	-L chrM:576-16024
 ```
 
 ### Keep variants located in the CR region
 
 ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk SelectVariants \
	-R Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	-V ${vcf_fiiltered_trimmed_shifted_chrM.simpleName}.vcf.gz \
	-O ${vcf_fiiltered_trimmed_shifted_chrM.simpleName}_ControlRegion.vcf.gz \
	-L chrM:7455-8576
  ```
  
  ### Shift the position of the variants located in the CR region back to the reference genome position
  
  ```
gzip -cd ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader.vcf.gz | awk ' \$2>=8001 {print \$1"\t"\$2-8000"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} \$2<=8000 {print \$1"\t"\$2+8569"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9"\t"\$10} ' | gzip -c  > ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_NoHeader_ShiftedBack.vcf.gz
```

### Sort the variants and index the file

```
bcftools sort ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1.vcf.gz -O z -o ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1_sorted.vcf.gz

singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk IndexFeatureFile \
        -I ${vcf_fiiltered_trimmed_CR_region_chrM.simpleName}_merged1_sorted.vcf.gz
```

## Merge the variants from the different samples and index file

```
bcftools merge -l ${list_MT_vcf} -O z -o ${version}_merged_chrM_Mutect2_filtered.vcf.gz

singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
        gatk --java-options "-Xmx4G" \
       IndexFeatureFile \
        -I ${version}_merged_chrM_Mutect2_filtered.vcf.gz
```

## Post variant calling QC

Cf main Readme file

## Frequency calculation for the MT variants
	
Only works for the MT frequency as it is necessary to calculate the VAF (Variant allele fraction or heteroplasmy levels)
	
```
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
```	

## MT variants annotation

Cf main Readme file
