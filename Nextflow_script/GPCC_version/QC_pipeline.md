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
  
### Mosdept   [Paper](https://pubmed.ncbi.nlm.nih.gov/29096012/)   -    [GitHub](https://github.com/brentp/mosdepth)
  

Command-line tool for rapidly calculating genome-wide sequencing coverage.
  

  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/mosdepth-0.3.2.sif \
mosdepth ${bam.simpleName} ${bam}
  ```
  
### Picard CollectWgsMetrics

Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
  
  ```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk	CollectWgsMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_collect_wgs_metrics.txt \
	-R ${ref_genome_cvmfs_file}
  ```

### PICARD CollectAlignmentSummaryMetrics

Produces a summary of alignment metrics from a SAM or BAM file. This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters.
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
gatk CollectAlignmentSummaryMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_Picard_Alignment
  ```

### PICARD QualityScoreDistribution

This tool is used for determining the overall 'quality' for a library in a given run. To that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases corresponding to those scores. 
  
  This one require R (Hence done while loading compute canada)
  
  ```
java -jar \$EBROOTPICARD/picard.jar \
	QualityScoreDistribution \
	I=${bam} \
	O=${bam.simpleName}_qual_score_dist.txt \
	CHART= ${bam.simpleName}_qual_score_dist.pdf
  ```

  PICARD full list : https://broadinstitute.github.io/picard/command-line-overview.html
  

### Bamtools Stats

**----- Bamtools Stats is not in the pipeline -----**

The command bamtools stats prints general alignment statistics from the BAM file.
  
  https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/data_manipulation_tools/bamtools/running_bamtools_commands/


## Agregation of Individual QC results using MultiQC

[MultiQC](https://multiqc.info) : Aggregate results from bioinformatics analyses across many samples into a single report
 
 Require to do one for the individual results (pre and post alignement) and another one for the population data

  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
multiqc ${params.outdir_ind}/${version}/QC/
  ```

## Vcf QC (Post-alignment)

### BcfTools stats

Parses VCF or BCF and produces text file stats which is suitable for machine processing and can be plotted using plot-vcfstats.
  
  https://samtools.github.io/bcftools/bcftools.html#stats
  
  ```
  bcftools stats  ${vcf_file}
  ```

### VCFTools TsTv-by-count

Calculates the Transition / Transversion ratio as a function of alternative allele count. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.count".

```
vcftools --gzvcf ${vcf_file} --TsTv-by-count --out ${vcf_file.simpleName}_Vcftools_TsTv_count
```

### VcfTools TsTv-by-qual

Calculates the Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.qual".

```
vcftools --gzvcf ${vcf_file} --TsTv-by-qual --out ${vcf_file.simpleName}_Vcftools_TsTv_qual
```
  
  VcfTools details : https://vcftools.github.io/man_latest.html

### VEP  statistics

VEP writes an HTML file containing statistics pertaining to the results of your job - Included by default, just redirected the file to the QC folder 
  
  https://m.ensembl.org/info/docs/tools/vep/vep_formats.html#stats
  
```
vep \
        -i ${vcf_file} \
        -o ${vcf_file.simpleName}_${version}_annotation_tab.tsv \
        --cache \
        --dir_cache /mnt/common/DATABASES/REFERENCES/GRCh38/VEP/ \
        --everything \
        --tab \
        --stats_file ${vcf_file.simpleName}_VEP_stats
```

## Agregator of population QC results using MultiQC
  
[MultiQC](https://multiqc.info) : Aggregate results from bioinformatics analyses across many samples into a single report
 
 Require to do one for the individual results (pre and post alignement) and another one for the population data


  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
multiqc $params.outdir_pop/${version}/QC/
```
  

