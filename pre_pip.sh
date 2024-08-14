#!/bin/bash

#SBATCH --job-name=pre_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --error=pre.%J.err     # standard error file
#SBATCH --output=pre.%J.out    # standard output file

#working directory
cd /g100_work/Gen4A_mulmy/
echo $PWD
sleep 5

#loading softwares
module load profile/bioinf
module load python

source snakemake_env/bin/activate

echo "pre-pipeline started"

snakemake --snakefile snakefiles/pre_hg19.snakefile --use-singularity

echo "The end"
