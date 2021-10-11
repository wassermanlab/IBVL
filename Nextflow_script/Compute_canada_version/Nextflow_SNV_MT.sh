#!/bin/sh

## CPU Usage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=200:30:00

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

module load singularity/3.7
module load nextflow
module load picard
module load bwa
module load samtools
module load bcftools

nextflow run Nextflow_SNV_MT_211009.nf -resume -with-trace
