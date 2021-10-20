
## Readme Goal

This readme is for everyone to comment on things to add / modify in the script


## .sh file

### Which tracing and vizualisation to use

Details are available here : https://www.nextflow.io/docs/latest/tracing.html

nextflow log <run name> --resume -with-report [file name] --with-trace -with-timeline [file name] -with-dag flowchart.png
  

## .config file

Full list available here : https://www.nextflow.io/docs/latest/config.html
  
process.executor = 'slurm'
process.queue = 'silent_q'

#Process.shell is necessary to avoid an error while unmounting GPCC and mounting CC modules through Nextflow
  
process.shell = ['/bin/bash','-e']

  
Add the limit in number of jobs submitted at a given time : 20?
queueSize : The number of tasks the executor will handle in a parallel manner (default: 100).

executor {
    queueSize = 20
}
  
## .nf file

### QC
  
  #### Fastq QC
  
  fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
  Would be great to use a second one?
  
  #### Bam QC

Picard CollectWgsMetrics      

Picard BamIndexStats  

PICARD CollectAlignmentSummaryMetrics   

PICARD QualityScoreDistribution   

mosdepth? 

#### SNV calls QC
  VEP ?

#### MT calls QC
  VEP?

#### SV callsQC
  ?

#### Agregator of QC results
  MultiQC
  
  ### Alignment (Fastq --> Bam)
  
  	bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${ref_genome_cvmfs_file} ${reads} | samtools view -Sb | samtools sort -o ${sampleId}_sorted.bam

  Piping the sam into bam and sorting allows to not save the sam file
  
  
  Bwa additional options : http://bio-bwa.sourceforge.net/bwa.shtml
  
  -t INT	Number of threads : Should be changed for GPCC?
  
  -R STR	Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. 
  
  
  Samtools options : http://www.htslib.org/doc/samtools-view.html
  
  -S The -S indicates the input is in SAM format and the "b" indicates that you'd like BAM output.
  Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.
  --> Try to remove it and compare outputs.
  
  
  ### SNV calling
  
  #### DeepVariant
  
  Deepvariant performs better than other tools such as GATK
  
  DeepVariant flags : https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
  
   #### GLnexus
  
  GLnexus is advised for joint calling following calling with DeepVariant
  
  GLnexus options : To find
  
  ### MT variant calling
  
  Decided to use GATK and follow the broad guidelines : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  It requires many step as reads mapping to the MT genome are extracted, and then re-aligned to a MT genome and a shifted MT genome to adress the circularity of the genome (And the reads mapping on the artificial breakpoint in the linear reference genome).
  
  The steps are described here : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-
  
  The general idea was kept while some step are slightly different
  
  ### Annotation of SNV and MT
  
  VEP options : https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html
  
  --everything = Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --var_synonyms, --variant_class, --mane
  
  --stats_file [filename] = Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"
