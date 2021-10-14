
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
  Bam : 

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
  
  
