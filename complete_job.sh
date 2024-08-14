#!/bin/bash

# Submit the pre-processing pipeline
job1_id=$(sbatch pre_pip.sh)
job1_id=$(echo $job1_id | awk '{print $4}')

# Submit the post-processing pipeline with dependency on the first
sbatch --dependency=afterok:$job1_id post_pip.sh

