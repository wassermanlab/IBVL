 # QC_pipeline file

This readme is for everyone to comment on things to add / modify in the script

This Readme only concerns the Quality control parts of the pipeline, including both pre and post alignement and post-calling QC.
  
  ## Fastq QC (pre-alignment)
  
### [fastqc]( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  
```
  fastqc ${fastq_file}
  ```
  
  ### Other tool
  
  Would be great to use a second one?
  
  ## Bam QC (Post-alignment)
  
### Mosdept

[Paper](https://pubmed.ncbi.nlm.nih.gov/29096012/)   -    [GitHub](https://github.com/brentp/mosdepth)

Command-line tool for rapidly calculating genome-wide sequencing coverage. 

  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/mosdepth-0.3.2.sif \
mosdepth ${bam.simpleName} ${bam}
  ```

### [PICARD](https://broadinstitute.github.io/picard/command-line-overview.html)

Picard offers several independant QC tools, the ones included in the pipeline are described below.

**[Picard CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-)**

Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
  
  ```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk	CollectWgsMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_collect_wgs_metrics.txt \
	-R ${ref_genome_cvmfs_file}
  ```

**[PICARD CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036883111-CollectAlignmentSummaryMetrics-Picard-)**

Produces a summary of alignment metrics from a SAM or BAM file. This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters.
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
gatk CollectAlignmentSummaryMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_Picard_Alignment
  ```

**[PICARD QualityScoreDistribution](https://gatk.broadinstitute.org/hc/en-us/articles/360037057312-QualityScoreDistribution-Picard-)**

This tool is used for determining the overall 'quality' for a library in a given run. To that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases corresponding to those scores. 
  
  This one require R (Hence done while loading compute canada)
  
  ```
java -jar \$EBROOTPICARD/picard.jar \
	QualityScoreDistribution \
	I=${bam} \
	O=${bam.simpleName}_qual_score_dist.txt \
	CHART= ${bam.simpleName}_qual_score_dist.pdf
  ```

### [Bamtools Stats](https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/data_manipulation_tools/bamtools/running_bamtools_commands/)

The command bamtools stats prints general alignment statistics from the BAM file.

**----- Bamtools Stats is not in the pipeline -----**
  
## Agregation of Individual QC results using MultiQC

[MultiQC](https://multiqc.info) : Aggregate results from bioinformatics analyses across many samples into a single report
 
 Require to do one for the individual results (pre and post alignement) and another one for the population data

  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
multiqc ${params.outdir_ind}/${version}/QC/
  ```

## Vcf QC (Post-alignment)

### [BcfTools stats](https://samtools.github.io/bcftools/bcftools.html#stats)

Used for SNV and MT vcf

Parses VCF or BCF and produces text file stats which is suitable for machine processing and can be plotted using plot-vcfstats.
  
  ```
bcftools stats ${vcf_file_MT/SNV} > ${vcf_file_MT/SNV.simpleName}_bcftools_stat
```

### [VcfTools](https://vcftools.github.io/man_latest.html)

VcfTools offers several QC options, the ones included in the pipeline are described below

**VcfTools TsTv-by-qual**

Calculates the Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.qual".

Used for SNV and MT vcf

```
vcftools --gzvcf ${vcf_file_MT/SNV} --TsTv-by-qual --out ${vcf_file_MT/SNV.simpleName}_Vcftools_TsTv_qual
```

**VCFTools TsTv-by-count**

Calculates the Transition / Transversion ratio as a function of alternative allele count. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.count".

**----- VCFTools TsTv-by-count is not in the pipeline -----**

Do not work on MT vcf : Error: Polyploidy found, and not supported by vcftools for MT variants

For SNV vcf, it does work on Genome Canada but not in GPCC (Even whle  unmounting GPCC and mounting Compute Canada)

Error message : Segmentation fault (core dumped)

### [VEP  statistics](https://m.ensembl.org/info/docs/tools/vep/vep_formats.html#stats)

VEP writes an HTML file containing statistics pertaining to the results of your job - Included by default, just redirected the file to the QC folder 

Used for SNV and MT vcf
  
```
vep \
        -i ${vcf_file_MT/SNV} \
        -o ${vcf_file_MT/SNV.simpleName}_${version}_annotation_tab.tsv \
	--offline \
        --cache \
        --dir_cache /mnt/common/DATABASES/REFERENCES/GRCh38/VEP/ \
        --everything \
        --tab \
        --stats_file ${vcf_file_MT/SNV.simpleName}_VEP_stats

```

## Agregator of population QC results using MultiQC
  
[MultiQC](https://multiqc.info) : Aggregate results from bioinformatics analyses across many samples into a single report
 
 Require to do one for the individual results (pre and post alignement) and another one for the population data


  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
multiqc $params.outdir_pop/${version}/QC/
```
  

