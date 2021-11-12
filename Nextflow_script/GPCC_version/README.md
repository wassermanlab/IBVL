
# Readme file

This readme is for everyone to comment on things to add / modify in the script

The complete pipeline (form fastq to annotated variant frequencies) is implemented within one nextflow pipeline.

However, for clarity and easier improvements, the pipeline description is splitted in several files :

- The [Alignment_SNV_pipeline.md](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/Alignment_SNV_pipeline.md) for the alignement step, SNV calling, SNV annotation and SNV frequency calculation. Some parts are common with the Mitochondrial and SV pipeline.

- The [Mitochondria_pipeline.md](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/Mitochondria_pipeline.md) for the part that are specific to the mitochondrial pipeline.

- The [QC_pipeline.md](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/QC_pipeline.md) for the quality control steps of the pipeline, including pre and post alignemnt QC.

- The SV_pipeline.md : To do.

# Nextflow config and sh files

For nextflow to work, a config file must be created, and the nextflow script must be launched using slurm scheduler and a sh file.

## .sh file

Which tracing and vizualisation to use, details are available here : https://www.nextflow.io/docs/latest/tracing.html

The <assembly_version> should be GRCh37 or GRCh38. It will automatically select the appropriate reference genome, VEP cache and mitochondrial nomenclature. It is a mendatory information.

```
nextflow log <run name> --resume -profile <assembly_version> -with-report [file name] --with-trace -with-timeline [file name] -with-dag flowchart.png
```  

## .config file

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

The profiles allows to guide the script to the appropriate reference genome, VEP cache and mitochondrial nomenclature. GRCh37 or GRCh38 must be specified in the .sh file

```
profiles {

        GRCh37 {
		params.assembly='GRCh37'
                params.ref_genome='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa'
                params.ref_genome_fai='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.fai'
                params.ref_genome_amb='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.amb'
                params.ref_genome_ann='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.ann'
                params.ref_genome_bwt='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.bwt'
                params.ref_genome_pac='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.pac'
                params.ref_genome_sa='/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa.sa'            
                params.Mitochondrial_chromosome='MT'
		params.vep_cache='/mnt/common/DATABASES/REFERENCES/GRCh37/VEP/'
        }

        GRCh38 {
		params.assembly='GRCh38'
                params.ref_genome='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa'
                params.ref_genome_fai='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai'
                params.ref_genome_amb='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb'
                params.ref_genome_ann='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann'
                params.ref_genome_bwt='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt'
                params.ref_genome_pac='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac'
                params.ref_genome_sa='/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa' 
		params.Mitochondrial_chromosome='chrM'
		params.vep_cache='/mnt/common/DATABASES/REFERENCES/GRCh38/VEP/'
        }
}
```



# .nf file

The .nf file contains the whole pipeline used to produce the IBVL. Each part of the pipeline is described in the other md files.
 
 # Pipeline overview
 
 Here is a simplified overview of the pipeline :
 
 
 <img width="922" alt="image" src="https://user-images.githubusercontent.com/54953390/141167242-8bef9c97-ab39-403f-b758-19708b811c12.png">

