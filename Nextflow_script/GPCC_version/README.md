
# Readme Goal

This readme is for everyone to comment on things to add / modify in the script


# .sh file

Which tracing and vizualisation to use, details are available here : https://www.nextflow.io/docs/latest/tracing.html

```
nextflow log <run name> --resume -with-report [file name] --with-trace -with-timeline [file name] -with-dag flowchart.png
```  

# .config file

Full list available here : https://www.nextflow.io/docs/latest/config.html

```
process.executor = 'slurm'
process.queue = 'silent_q'
```

#Process.shell is necessary to avoid an error while unmounting GPCC and mounting CC modules through Nextflow
```  
process.shell = ['/bin/bash','-e']
```
  
Add the limit in number of jobs submitted at a given time : 20?
queueSize : The number of tasks the executor will handle in a parallel manner (default: 100).
```
executor {
    queueSize = 20
}
```  
# .nf file

  
  ## Alignment (Fastq --> Bam)
  
```
bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${ref_genome_cvmfs_file} ${reads} | samtools view -Sb | samtools sort -o ${sampleId}_sorted.bam`
samtools index ${sampleId}_sorted.bam
```
  
  
  Piping the sam into bam and sorting allows to not save the sam file
  
  
  Bwa additional options : http://bio-bwa.sourceforge.net/bwa.shtml
  
  -t INT	Number of threads : Should be changed for GPCC?
  
  -R STR	Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. 
  
  
  Samtools options : http://www.htslib.org/doc/samtools-view.html
  
  -S The -S indicates the input is in SAM format and the "b" indicates that you'd like BAM output.
  Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.
  --> Try to remove it and compare outputs.
  
  
  ## SNV calling
  
  ### DeepVariant
  
  Deepvariant performs better than other tools such as GATK
  
  DeepVariant flags : https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
	
```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ -B /mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/ /mnt/common/SILENT/Act3/singularity/deepvariant-1.2.0.sif \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${ref_genome_file} \
	--reads=${bam.simpleName}.bam \
	--regions chr20 \
	--output_gvcf=${bam.simpleName}.g.vcf.gz \
	--output_vcf=${bam.simpleName}.vcf.gz
```
  
   ### GLnexus
  
  GLnexus is advised for joint calling following calling with DeepVariant
  
  GLnexus options : To find
	
```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/glnexus-1.4.1.sif \
	glnexus_cli \
	--config DeepVariant \
	--list ${list_gvcf} > DeepVariant_GLnexus_${version}.bcf
	
bcftools view ${bcf_file} | bgzip -c > ${bcf_file.simpleName}.vcf.gz

singularity exec -B /home -B /project -B /scratch -B /localscratch /home/correard/scratch/SnakeMake_VariantCalling/MT/Mutect2/gatk.sif \
        gatk --java-options "-Xmx4G" \
       IndexFeatureFile \
        -I ${bcf_file.simpleName}.vcf.gz
```
  
  ## MT variant calling
  
  Decided to use GATK and follow the broad guidelines : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  It requires many step as reads mapping to the MT genome are extracted, and then re-aligned to a MT genome and a shifted MT genome to adress the circularity of the genome (And the reads mapping on the artificial breakpoint in the linear reference genome).
  
  The steps are described here : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  The general idea was kept while some step are slightly different
  
  ## Frequency calculation and annotation of SNV and MT variants
  
###Frequency calculation for the SNV 
	
Only works for the SNV frequency as other values are necessary for the MT variants
	
```
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
```
	
### Frequency calculation for the MT variants
	
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
	
	### Annotation for SNV and MT variants
	
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
  
  VEP options : https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html
  
  --everything = Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --var_synonyms, --variant_class, --mane
  
  --stats_file [filename] = Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"

  
  ## QC
  
  ### Fastq QC (pre-alignment)
  
  **fastqc** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
```
  fastqc ${fastq_file}
  ```
  
  Would be great to use a second one?
  
  ### Bam QC (Post-alignment)
  
**Picard CollectWgsMetrics** : Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
  
  ```
singularity exec -B /mnt/common/DATABASES/REFERENCES/ -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk	CollectWgsMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_collect_wgs_metrics.txt \
	-R ${ref_genome_cvmfs_file}
  ```

**PICARD CollectAlignmentSummaryMetrics**   : Produces a summary of alignment metrics from a SAM or BAM file. This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters.
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
gatk CollectAlignmentSummaryMetrics \
	-I ${bam} \
	-O ${bam.simpleName}_Picard_Alignment
  ```

**PICARD QualityScoreDistribution** : This tool is used for determining the overall 'quality' for a library in a given run. To that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases corresponding to those scores. 
  
  This one require R (Hence done while loading compute canada)
  
  ```
java -jar \$EBROOTPICARD/picard.jar \
	QualityScoreDistribution \
	I=${bam} \
	O=${bam.simpleName}_qual_score_dist.txt \
	CHART= ${bam.simpleName}_qual_score_dist.pdf
  ```

  PICARD full list : https://broadinstitute.github.io/picard/command-line-overview.html
  
**mosdepth** : command-line tool for rapidly calculating genome-wide sequencing coverage.
  
  Paper : https://pubmed.ncbi.nlm.nih.gov/29096012/
  GitHub : https://github.com/brentp/mosdepth
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/mosdepth-0.3.2.sif \
mosdepth ${bam.simpleName} ${bam}
  ```

**Bamtools Stats** : The command bamtools stats prints general alignment statistics from the BAM file.
  
  https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/data_manipulation_tools/bamtools/running_bamtools_commands/

### Vcf QC (Post-alignment)[SNV, MT, SV]

**BcfTools stats**   : Parses VCF or BCF and produces text file stats which is suitable for machine processing and can be plotted using plot-vcfstats.
  
  https://samtools.github.io/bcftools/bcftools.html#stats
  
  ```
  bcftools stats  ${vcf_file}
  ```

**VCFTools TsTv-by-count**   : Calculates the Transition / Transversion ratio as a function of alternative allele count. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.count".


**VcfTools TsTv-by-qual**  : Calculates the Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.qual".

```
vcftools --vcf ${vcf_file} --TsTv-by-count --TsTv-by-qual --out ${vcf_file}_Bcftools_stats
```
  
  VcfTools details : https://vcftools.github.io/man_latest.html

**VEP  statistics** : VEP writes an HTML file containing statistics pertaining to the results of your job - Included by default, just redirected the file to the QC folder 
  
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

### Agregator of QC results
  
 **MultiQC** : Aggregate results from bioinformatics analyses across many samples into a single report
  
  ```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/multiqc-1.9.sif \
multiqc ${params.outdir_ind}/${version}/QC/
  ```
  
  https://multiqc.info
