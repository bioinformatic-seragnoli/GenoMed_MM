configfile: "config_files/config_hg38.yaml"
configfile: "config_files/samples_list.yaml"

rule all:
    input:
        expand("results/BOBaFIT/{tumor}_chrarmClass.tsv", tumor=config["bam"]),
        expand("results/BOBaFIT/{tumor}_BOB_seg.tsv", tumor=config["bam"]),
        expand("results/ichorCNA/{tumor}.cna.seg", tumor=config["bam"]),
        expand("results/readDepth/{tumor}.bin{binSize}.wig", tumor=config["bam"], binSize=str(config["binSize"]))

rule read_counter:
    input:
        lambda wildcards: config["bam"][wildcards.tumor]
    output:
        temp("results/readDepth/{tumor}.bin{binSize}.wig")
    params:
        wd=config["wd"],
        readCounter=config["readCounterScript"],
        binSize=config["binSize"],
        qual="20",
        chrs=config["chrs"]
    log:
        "logs/readDepth/{tumor}.bin{binSize}.log"
    shell:
        "singularity run -B {params.wd} singularity/hmm_utils_1.0.0.sif {params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
    input:
        tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig"
    output:
        cna="results/ichorCNA/{tumor}.cna.seg",
        segTxt="results/ichorCNA/{tumor}.seg.txt"
    threads: 5
    params:
        wd=config["wd"],
        outDir="results/ichorCNA/",
        rscript=config["ichorCNA_rscript"],
        id="{tumor}",
        ploidy=config["ichorCNA_ploidy"],
        normal=config["ichorCNA_normal"],
        gcwig=config["ichorCNA_gcWig"],
        mapwig=config["ichorCNA_mapWig"],
        normalpanel=config["ichorCNA_normalPanel"],
        estimateNormal=config["ichorCNA_estimateNormal"],
        estimatePloidy=config["ichorCNA_estimatePloidy"],
        estimateClonality=config["ichorCNA_estimateClonality"],
        scStates=config["ichorCNA_scStates"],
        maxCN=config["ichorCNA_maxCN"],
        includeHOMD=config["ichorCNA_includeHOMD"],
        chrs=config["ichorCNA_chrs"],
        chrTrain=config["ichorCNA_chrTrain"],
        genomeBuild=config["ichorCNA_genomeBuild"],
        genomeStyle=config["ichorCNA_genomeStyle"],
        centromere=config["ichorCNA_centromere"],
        fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=config["ichorCNA_minMapScore"],
        maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
        exons=config["ichorCNA_exons"],
        txnE=config["ichorCNA_txnE"],
        txnStrength=config["ichorCNA_txnStrength"],
        plotFileType=config["ichorCNA_plotFileType"],
        plotYlim=config["ichorCNA_plotYlim"],
        libdir=config["ichorCNA_libdir"]
    log:
        "logs/ichorCNA/{tumor}.log"    
    shell:
        "singularity run -B {params.wd} singularity/ichorcna_1.0.0.sif {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

rule BOBaFIT:
    input:
        seg="results/ichorCNA/{tumor}.seg.txt"
    output:
        bobseg="results/BOBaFIT/{tumor}_BOB_seg.tsv",
        bincall="results/BOBaFIT/{tumor}_chrarmClass.tsv"
    params:
        wd=config["wd"],
        outDir="results/BOBaFIT/",
        rscript=config["BOBaFIT_rscript"],
        normChr=config["BOBaFIT_normChr"],
        refGenome=config["BOBaFIT_refGenome"],
        libdir=config["BOBaFIT_libdir"]
    log:
        "logs/BOBaFIT/{tumor}.log"
    shell:
        "singularity run -B {params.wd} singularity/bobafit_1.0.0.sif {params.rscript} --input {input.seg} --normChr \"{params.normChr}\" --refGenome {params.refGenome} --libdir {params.libdir} --output {params.outDir} > {log} 2> {log}"

		
