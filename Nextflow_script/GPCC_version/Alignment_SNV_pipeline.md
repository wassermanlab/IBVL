 # Alignment_SNV_pipeline file

This readme is for everyone to comment on things to add / modify in the script

This Readme concerns the alignement and SNV part of the pipeline. Some parts are common with the MT and the SV pipeline.
 
 **GRCh37 and GRCh38**
 
 The first release of the IBVL will be done based on GRCh37, however, it is the team wishes to also made a release available in GRCh38 quickly after.
 
 Changing the assembly requires to change the reference genome for alignment, as well as the VEP cache for annotation.
 
 The user should comment out with '\\' the version they do not want the pipeline to use
 
 ```
// For GRCh38 - Get the genome files
// To update : When using GRCh38, should use "chrM" to extract MT reads and chr20 to extract only variants on chr 20 with DeepVariant 
//params.ref_genome="/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa"
//params.vep_cache="/mnt/common/DATABASES/REFERENCES/GRCh38/VEP/"

// For GRCh37 - Get the genome files
// To update : When using GRCh38, should use "MT" to extract MT reads and 20 to extract only variants on chr 20 with DeepVariant
params.ref_genome="/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa"
params.vep_cache="/mnt/common/DATABASES/REFERENCES/GRCh37/VEP/"

ref_genome_file = file(params.ref_genome)
vep_cache=file(params.vep_cache)
```
 
 **Future IBVL update:** The user should define which assembly they want to use when lauching the workflow (nextflow run <pipeline-name> --genome GRCh37 <--other-arguments>
	
Need to modify the config file with something in the lines of :
	
 ```
params {
  genomes {
    'GRCh37' {
		ref_genome="/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa"
		vep_cache="/mnt/common/DATABASES/REFERENCES/GRCh37/VEP/"
		Mitochondiral_chromosome="chrM"
    }
    'GRCh38' {
		ref_genome="/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa"
		vep_cache="/mnt/common/DATABASES/REFERENCES/GRCh38/VEP/"
		Mitochondiral_chromosome="MT"
	}
}
 ```
 
 ## Alignment (Fastq --> Bam)
  
```
bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${ref_genome_file} ${reads} | samtools view -Sb | samtools sort -o ${sampleId}_sorted.bam`
samtools index ${sampleId}_sorted.bam
```
  
  Piping the sam into bam and sorting allows to not save the sam file
    
  **[Bwa](http://bio-bwa.sourceforge.net/bwa.shtml)**
  
  -t INT	Number of threads : Should be changed for GPCC?
  
  -R STR	Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. 
  
  
  **[Samtools](http://www.htslib.org/doc/samtools-view.html)**
  
  -S The -S indicates the input is in SAM format and the "b" indicates that you'd like BAM output.
  Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.
  
  ## SNV calling
  
  ### [DeepVariant](https://github.com/google/deepvariant)
  
  Deepvariant performs better than other tools such as GATK
  
  [DeepVariant flags](https://cloud.google.com/life-sciences/docs/tutorials/deepvariant)
	
```
singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ -B /mnt/common/DATABASES/REFERENCES/ /mnt/common/SILENT/Act3/singularity/deepvariant-1.2.0.sif \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${ref_genome_file} \
	--reads=${bam.simpleName}.bam \
	--regions chr20 \
	--output_gvcf=${bam.simpleName}.g.vcf.gz \
	--output_vcf=${bam.simpleName}.vcf.gz
```
  
  **Remove the --regions chr20 tag when runninig the final dry run**
  
   ### [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
  
  GLnexus is advised for joint calling following calling with DeepVariant
	
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
  
  
  ## Frequency calculation for the SNV using [VariantsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable)
	
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
        -F MULTI-ALLELIC 
```
	

	
## Annotation for SNV and MT variants using [VEP](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)

Common for SNV and MT variants
	
```
	vep \
        -i ${vcf_file_MT/SNV} \
        -o ${vcf_file_MT/SNV.simpleName}_${version}_annotation_tab.tsv \
	--offline \
        --cache \
        --dir_cache ${vep_cache} \
        --everything \
        --tab \
        --stats_file ${vcf_file_MT/SNV.simpleName}_VEP_stats
```
  
  --everything = Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --var_synonyms, --variant_class, --mane
  
  --stats_file [filename] = Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"

**Future IBVL Update :** : Include VEP plugins
