configfile: "config_files/config_hg19.yaml"
configfile: "config_files/samples_list.yaml"

rule all:
    input:
        expand("results/remDup_recal_bam/remDup_recal_{tumor}.bam.bai", tumor=config["samples"]),
        expand("results/remDup_recal_bam/remDup_recal_{tumor}.bam", tumor=config["samples"]),
        expand("results/removeDups/remDup_{tumor}.bam", tumor=config["samples"]),
        expand("results/alignedFastq/{tumor}.bam", tumor=config["samples"]),
        expand("results/alignedFastq/{tumor}.sam", tumor=config["samples"]),
        expand("results/mergedFastq/{samples}.fastq.gz", samples=config["samples"])

rule merge_fastq:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        "results/mergedFastq/{samples}.fastq.gz"
    log:
        "logs/mergedFastq/{samples}.log"
    shell:
        "cat {input} > {output} 2> {log}"

rule align_fastq:
    input:
        reads="results/mergedFastq/{tumor}.fastq.gz"
    output:
        temp("results/alignedFastq/{tumor}.sam")
    params:
        wd=config["wd"],
        rg=r"@RG\tID:" + str(config["bwa_id"]) + r"\tLB:{tumor}\tSM:{tumor}\tPL:illumina",
        refGen=config["refGen"]
    threads: 4
    benchmark:
        "results/benchmarks/{tumor}_bwa_benchmark.txt"
    log:
        "logs/alignedFastq/{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/bwa_1.0.0.sif mem -R '{params.rg}' -t {threads} {params.refGen} {input.reads} -o {output} 2> {log}"

rule sort_sam:
    input:
        sam="results/alignedFastq/{tumor}.sam"
    output:
        temp("results/alignedFastq/{tumor}.bam")
    params:
        wd=config["wd"]
    log:
        "logs/alignedFastq/sort_{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/samtools_1.0.0.sif sort {input} -O bam -o {output} 2> {log}"

rule mark_dups:
    input:
        bam="results/alignedFastq/{tumor}.bam"
    output:
        temp("results/removeDups/remDup_{tumor}.bam")
    params:
        wd=config["wd"]
    log:
        mx="results/metrics/removeDups/{tumor}_remDupMetrics.txt",
        mark="logs/removeDups/remove_duplicates_{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/picard_1.0.0.sif MarkDuplicates I={input.bam} O={output} M={log.mx} 2> {log.mark}"

rule index_bam:
    input:
        bam="results/removeDups/remDup_{tumor}.bam"
    output:
        "results/removeDups/remDup_{tumor}.bam.bai"
    params:
        wd=config["wd"]
    benchmark:
        "results/benchmarks/{tumor}_indexing_benchmark.txt"
    threads: 2
    log:
        "logs/indexing/index_{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/samtools_1.0.0.sif index {input.bam} 2> {log}"

rule realign_target_creator:
    input:
        bam="results/removeDups/remDup_{tumor}.bam"
    output:
        "results/realign_bam/intervals/remDup_{tumor}.intervals"
    params:
        wd=config["wd"],
        refGen=config["refGen"],
        indels=config["gold_indels"]
    benchmark:
        "results/benchmarks/{tumor}_localalign_benchmark.txt"
    log:
        "logs/gatk/realign_target_creator_{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/gatk3_1.0.0.sif RealignerTargetCreator -R {params.refGen} -I {input.bam} --known {params.indels} -o {output} 2> {log}"

rule base_recalibrator:
    input:
        bam="results/removeDups/remDup_{tumor}.bam"
    output:
        "results/metrics/recal_tables/{tumor}_recalTable.table"
    log:
        "logs/gatk/base_recalibrator_{tumor}.log"
    params:
        wd=config["wd"],
        refGen=config["refGen"],
        indels=config["gold_indels"],
        gnomad=config["gnomad_af"]
    threads: 4
    benchmark:
        "results/benchmarks/{tumor}_baserecalibrator_benchmark.txt"
    shell:
        "singularity run -B {params.wd} singularity/gatk_4.3.0.0.sif gatk BaseRecalibrator -R {params.refGen} -I {input.bam} --known-sites {params.indels} --known-sites {params.gnomad} -O {output} 2> {log}"

rule apply_BQSR:
    input:
        bam="results/removeDups/remDup_{tumor}.bam",
        table="results/metrics/recal_tables/{tumor}_recalTable.table"
    output:
        "results/remDup_recal_bam/remDup_recal_{tumor}.bam"
    log:
        "logs/gatk/apply_BQSR_{tumor}.log"
    params:
        wd=config["wd"],
        refGen=config["refGen"]
    threads: 4
    benchmark:
        "results/benchmarks/{tumor}_applybqsr_benchmark.txt"
    shell:
        "singularity run -B {params.wd} singularity/gatk_4.3.0.0.sif gatk ApplyBQSR -R {params.refGen} -I {input.bam} -bqsr {input.table} -O {output} 2> {log}"

rule index_remdup_bam:
    input:
        bam="results/remDup_recal_bam/remDup_recal_{tumor}.bam"
    output:
        "results/remDup_recal_bam/remDup_recal_{tumor}.bam.bai"
    params:
        wd=config["wd"]
    threads: 2
    log:
        "logs/indexing/index_recal_{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/samtools_1.0.0.sif index {input.bam} 2> {log}"
