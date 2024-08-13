# GENOMED4ALL - GENOMIC PIPELINES STORAGE T5.2 | MULTIPLE MYELOMA PIPELINE 

## Aim of the pipeline
The following pipeline aims to classify copy number alterations (CNA) in Multiple Myeloma (MM) samples. The pipeline is optimized to classify these genomic alterations at chromosome arm level on 
ULP-WGS-sequenced genomic DNA (gDNA) samples. Specifically, we focus only on clonal CNAs, i.e. present in at least 50% of cancer cells.  

## Structure of the pipeline
The entire pipeline can be divided into two steps, a pre-processing and a post-processing one.

	* Pre-processing step: preparatory step of raw data for the downstream analyses (from fastq to bam files). Involved tools: samtools (sorting, indexing), bwa (alignment), picard (duplicates removal), 
	gatk (local realignment, base recalibration). Programming languages: bash, java. Computationally-expensive part of the pipeline. 
	* Post-processing step: customized part of the pipeline, which aim to categorize CN calls (from bam files to final output files). Involved tools: readCounter, ichorCNA, BOBaFIT. Languages: bash, R. 

## How to run a job on a single sample
All bioinformatic tools involved are run as singularity containers which are concatenated in a linear manner with Snakemake. Management of each job, i.e. analysis of a single sample per time, is 
orchestrated using SLURM. 
To adapt the pipeline to each user's directories refer to config files present in the config_files/ folder. samples_list.yaml must contain the correct path to the input sample files (either 
fastq or bam files); config_hg19.yaml and config_hg38.yaml contain both parameters necessary to correctly run all tools and thus can be modified depending on each user's needs (strictly in bash). 

To run the complete pipeline:
	sbatch complete_job.sh #to submit the pipeline
	squeue -u <username> #checking job progressing status

Both human reference genomes, i.e. hg19, hg38, are supported:
	vim pre_pip.sh #open the corresponding bash script
	snakemake --snakefile snakefiles/pre_hg19.snakefile --use-singularity #change the snakefile name according to the reference you want to test

The same steps must be followed for post_pip.sh. 

There is also the possibility to exclusively run either the pre- or post-processing steps of the pipeline, assuming that the required input files are present in the specified folders. 
To run the pre-processing step:
 	sbatch pre_pip.sh #to run only the pre-processing step (required fastq files must be located in the input/ folder)

To run the post-processing step:
	sbatch post_pip.sh #to run only the post-processing step (required bam files must be located in the results/ folder)

## Output
All tools' output files are located inside the results/ folder, organized into subfolders allocating non-temporary files. The final output of the entire pipeline is the <sample_name>_chrarmClass.tsv file located inside the subfolder BOBaFIT/, which reports the weighted CN value and corresponding binary classification (amplified or deleted) of each chromosome arm. 

## References
Singularity containers were directly pulled from DockerHub and build from custom Dockerfiles.
For more info, go to the following GitHub repository: https://github.com/gaia-mazzocchetti/GenoMed_MM
(ora privato, mettere pubblico? sul seragnoli? toglierlo del tutto? boh) 





