
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

```
nextflow log <run name> --resume -with-report [file name] --with-trace -with-timeline [file name] -with-dag flowchart.png
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
# .nf file

The .nf file contains the whole pipeline used to produce the IBVL. Each part of the pipeline is described in the other md files.
 
 # Pipeline overview
 
 Here is a simplified overview of the pipeline :
 
 
 <img width="922" alt="image" src="https://user-images.githubusercontent.com/54953390/141167242-8bef9c97-ab39-403f-b758-19708b811c12.png">

