 # Alignment_SNV_pipeline file

This readme is for everyone to comment on things to add / modify in the script

This Readme concerns the alignement and SNV part of the pipeline. Some parts are common with the MT and the SV pipeline.
 
 
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
  
  
  ## Frequency calculation for the SNV 
	
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
	

	
## Annotation for SNV and MT variants

Common for SNV and MT variants
	
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

