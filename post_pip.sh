#!/bin/bash

#SBATCH --job-name=post_pipeline
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=00:40:00
#SBATCH --error=post.%J.err     # standard error file
#SBATCH --output=post.%J.out    # standard output file

#working directory
cd /g100_work/Gen4A_mulmy/
echo $PWD
sleep 5

#loading softwares
module load profile/bioinf
module load python

source snakemake_env/bin/activate

echo "post-pipeline started"

snakemake --snakefile snakefiles/post_hg19.snakefile --use-singularity 

echo "The end! yey"
