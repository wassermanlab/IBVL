
# Mitochondria_pipeline file

This readme is for everyone to comment on things to add / modify in the script

This Readme only concerns the Mitochondrial part of the pipeline. Some parts are common with the SNV pipeline and therefore are described in the [alignment_SNV_pipeline](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/Alignment_SNV_pipeline.md) document.

## References downloading and indexing using [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)

It is necessary to download the reference genome as well as the shifted reference genomes

Fasta downloaded from [here](https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references)

```
bwa index ${ref_genome_MT_file}
bwa index ${ref_genome_MT_shifted_file}
 ```

## MT variant calling
  
The mitochondrial genome poses several challenges to the identification and understanding of somatic variants. To address the different challenges, the IBVL MT pipeline was created using several sources of information, manly :

- [The GATK best practices workflow for Mitochondrial short variant discovery (SNVs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-). Most of the text in the repo is coming from the GATK website.

- [Laricchia et al, 2021, bioRXiv](https://www.biorxiv.org/content/10.1101/2021.07.23.453510v1.full.pdf)
  
The pipeline is composed of many step described below. The general idea from the GATK pipeline was kept while some step are slightly different as the full description from GATK was not avalable at the time of implementation of the IBVL pipeline.
 
 ### Subset to keep only the reads mapped to the mitochondrial genome using [PrintReads](https://gatk.broadinstitute.org/hc/en-us/articles/360036715871-PrintReads)
 
 From GATK best practices website : This step filters the input WGS file to keep only the ChrM mapped reads. Tools involved: PrintReads - SubsetBamtoChrM

From BioRXiv paper : Pull reads from chrM (GATK PrintReads --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter)

 **Potential future update in the IBVL pipeline** : Also extract the reads that maps to NuMTs.
 
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
  
  ### Map the chrM reads to the MT reference genome and the shifted MT reference genome
  
From GATK best practices website : Revert the ChrM mapped reads from an aligned BAM to an unaligned BAM file. Tools involved: RevertSam (This step reverts the aligned BAM file containing only the reads mapped to the mitochondria to remove all alignment information while retaining the recalibrated base qualities and original alignment tags.) and Align the unmapped BAM file with the reference aligned BAM and shifted reference aligned BAM. Tools involved: MergeBamAlignment (This step merges the unaligned BAM file with the reference BAM file for the mitochondrial genome. The BAM file must also be aligned with the shifted mitochondrial BAM file. The shifted reference moves the breakpoint of the mitochondrial genome from the non-coding control region to the opposite side of the contig. This allows for sensitivity in the control region to account for variability across individuals.)

From BioRXiv paper : Not mentionned
    
  For the IBVL pipeline, the ChrM mapped reads from an aligned BAM are converted back to fastq (losing the recalibrated base qualities and original alignment tags) and then mapped again the reference genome for the mitochondiral genome and the shifted version of the mitochondrial genome.
  

**Convert back the bam file to fastq using [SamToFastq](https://gatk.broadinstitute.org/hc/en-us/articles/360036485372-SamToFastq-Picard-)**

```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk SamToFastq \
	INPUT=${chr_bam.baseName}.bam \
	FASTQ=${chr_bam.baseName}.fastq \
	INTERLEAVE=true \
	NON_PF=true

 ```
 
**Align Fastq to MT reference genome (also run using the shifted MT reference genome) using [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) and [SamTools](http://www.htslib.org/doc/samtools-index.html)**
 
 ```
bwa mem -R "@RG\\tID:${fastqfromsam.baseName}\\tSM:${fastqfromsam.baseName}\\tPL:illumina" Homo_sapiens_assembly38.chrM.fasta ${fastqfromsam.baseName}.fastq | samtools view -u -bS | samtools sort > ${fastqfromsam.baseName}_chrM.bam
  
samtools index ${fastqfromsam.baseName}_chrM.bam
 ```
 
### Identify Duplicate Reads using [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360036729891-MarkDuplicates-Picard-) for both bam files mapped with the reference genome and the shifted reference genome

From GATK best practices website : Tools involved: MarkDuplicates. This step identifies and tags duplicate reads in the aligned BAM files.

From BioRXiv paper : Considered done for the nuclear genomes : Briefly, GATK version 4.1.2.0 (McKenna et al. 2010) tools were used to estimate the median nuclear genome coverage (Picard CollectWgsMetrics), to exclude duplicates (Picard MarkDuplicates)

```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk MarkDuplicates \
      I=${bam/shifted_MT} \
      O=${bam/shifted_MT.simpleName}_marked_duplicates.bam \
      M=${bam/shifted_MT.simpleName}_marked_duplicates_metrics.txt

samtools index ${bam_MT.simpleName}_marked_duplicates.bam
```

### Collect coverage and performance metrics for BAM file using [CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036804671-CollectWgsMetrics-Picard-) for both bam files mapped with the reference genome and the shifted reference genome

From GATK best practices website : Tools involved: CollectWgsMetrics. This step collects several metrics for evaluating the coverage and performance of the WGS experiment. This includes the fraction of bases that pass base and mapping quality filters as well as coverage levels.

From BioRXiv paper : Considered done for the nuclear genomes : Briefly, GATK version 4.1.2.0 (McKenna et al. 2010) tools were used to estimate the median nuclear genome coverage (Picard CollectWgsMetrics), to exclude duplicates (Picard MarkDuplicates)

```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk CollectWgsMetrics \
        -I ${bam/shifted_MT} \
        -O ${bam/shifted_MT.simpleName}_collect_wgs_metrics_MT.txt \
        -R ${ref_genome_MT/shifted_file}
```
 
 ### Call MT variants in aligned BAM files using [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360051306691-Mutect2) for both bam files mapped with the reference genome and the shifted reference genome
 
This step is performed using the MT reeference genome and the shifted MT reference genome.

From GATK best practices website : Tools involved: Mutect2 - This step calls mitochondrial variants in the non-control region using the BAM file aligned to the mitochondrial reference and in the control region using the BAM file aligned to the shifted mitochondrial reference. Running Mutect2 in mitochondrial mode automatically sets parameters appropriately for calling on mitochondria with the --mitochondria flag. Specifically, the mode sets –-initial-tumor-lod to 0, –-tumor-lod-to-emit to 0, --af-of-alleles-not-in-resource to 4e-3, and the advanced parameter --pruning-lod-threshold to -4. Setting the advanced option --median-autosomal-coverage argument (default 0) activates a recommended filter against likely erroneously mapped NuMTs (nuclear mitochondrial DNA segments). For the value, provide the median coverage expected in autosomal regions with coverage.


From BioRXiv paper : Call variants (GATK Mutect2 --mitochondria-mode -- annotation StrandBiasBySample --max-reads-per-alignment-start 75 --max-mnp-distance 0).

**Future update in the IBVL pipeline** : "A USER ERROR has occurred: median-autosomal-coverage is not a recognized option"
Solenne asked GATK dev team if this flag will be included in future release of GATK4 or if they should not be used anymore
 
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
  
  ### Shift back the variants in the control region and merge the variants for one sample
  
  From GATK best practices website : 3 steps.

1. Liftover the output VCF files. Tools involved: LiftoverVcf - This step returns the variant calls back to the standard numbering system with the original alignment (OA) tags.
2. Combine the variant calls from the control region with the non-control region. Tools involved: MergeVcf - This step merges the output VCF file for the control region (BAM aligned to shifted reference) with the VCF file for the non-control region into a single variant file.
3. Merge stats files for output VCFs. Tools Involved: Mutect2- MergeMutectStats - This step merges the stats file for the variant calls of the control region with the stats file for the variant calls of the non-control region.

From BioRXiv paper : Variants called on the shifted reference were mapped back to standard coordinates (Picard LiftOver) and combined with variants from the non-control region.
  
  **LiftOver the variants called with the shifted reference using [LiftOverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360036347792-LiftoverVcf-Picard-)**
  
  ShiftBack.chain downloaded from [here](https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references)
  
  ```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk LiftoverVcf \
	I=${MT_call_variants_shifted} \
	O=${MT_call_variants_shifted.simpleName}_lifted_over.vcf \
	CHAIN=${ShiftBack_chain_MT_file} \
	REJECT=${MT_call_variants_shifted.simpleName}_rejected_variants.vcf \
	R=${ref_genome_MT_file}
```
  
**Combine the variant calls from the control region with the non-control region using [MergeVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard-)**
  
  ```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk MergeVcfs \
          I=${MT_call_variants} \
          I=\${sample_name}_sorted_chrM_chrM_shifted_marked_duplicates_Mutect2_lifted_over.vcf \
          O=\${sample_name}_merged.vcf.gz
```


**Merge stats files for output VCFs using [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360051306691-Mutect2)**

This step is necessary for the following steps of variant filetring (The stat file is essential for FilterMutectCalls

```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk MergeMutectStats \
        -stats ${MT_call_variants_stat} \
        -stats \${sample_name}_sorted_chrM_chrM_shifted_marked_duplicates_Mutect2.vcf.gz.stats \
        -O \${sample_name}_merged.stats
```
 
 ### Filter merged Mutect2 calls
 
 From GATK best practices website : Tools involved: FilterMutectCalls - This step filters the output VCF files based on specific parameters, such as a minimum allele fraction, maximum alternate allele count, and estimate of contamination. The --autosomal-coverage parameter specifically filters out potential NuMTs. Specifying the --mitochondrial-mode parameter automatically sets the filters to the mitochondrial defaults.

From BioRXiv paper : Mutect2 variants were then filtered (GATK FilterMutectCalls --stats raw_vcf_stats --max-alt-allele-count 4 --mitochondria-mode --autosomal_coverage nDNA_MEDIAN_COV --min_allele_fraction 0.01) and multi-allelic sites were split into different variants (LeftAlignAndTrimVariants --split-multi-allelics --dont-trim-alleles --keep-original-ac)

**Variant filetrin using [FilterMutectCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360036715731-FilterMutectCalls)**
 
 ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk FilterMutectCalls \
	-V ${MT_MergeVcfs.simpleName}.vcf.gz \
	-R Homo_sapiens_assembly38.chrM.fasta \
	--stats ${MT_MergeVcfs.simpleName}.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	-O ${MT_MergeVcfs.simpleName}_filtered.vcf.gz
```

**Future update in the IBVL pipeline** : "A USER ERROR has occurred: autosomal-coverage is not a recognized option"
Solenne asked GATK dev team if this flag will be included in future release of GATK4 or if they should not be used anymore

**Variant trimming using [LeftAlignAndTrimVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037225872-LeftAlignAndTrimVariants)**

```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk LeftAlignAndTrimVariants \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${MT_Filter_Mutect_Calls.simpleName}.vcf.gz \
	-O ${MT_Filter_Mutect_Calls.simpleName}_trimmed.vcf.gz \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
```
 
 ### Filter out Blacklisted Sites
 
**Future update in the IBVL pipeline** : Add filter out Blacklisted Sites - Tools involved: VariantFiltration. This step filters out blacklisted sites containing unwanted artifacts.

blacklist_sites.hg38.chrM.bed and blacklist_sites.hg38.chrM.bed.idx downloaded from [here](https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references)

```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk VariantFiltration \
   -R Homo_sapiens_assembly38.chrM.fasta \
   -V ${MT_trimmed.simpleName}.vcf.gz \
   -O ${MT_trimmed.simpleName}_filtered.vcf.gz \
   --mask-name "GATK_artifact" \
   --mask ${blacklist_sites.hg38.chrM}
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

Cf main Alignement and SNV calling file

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
        -GF AF
```	

## MT variants annotation

Cf main Alignement and SNV calling file
