
# Readme file

This readme is for everyone to comment on things to add / modify in the script

This readme file concern the part that are common for all the pipelines (alignemnt, pre and post alignement QC) and the SNV calling.

For the other parts of the pipeline, refer to :

- The [Readme_MT](https://github.com/scorreard/IBVL/blob/main/Nextflow_script/GPCC_version/Readme_MT.md) file for the parts specific to MT variant calling 

- The Readme_SV file for the parts specific to SV calling


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

 
