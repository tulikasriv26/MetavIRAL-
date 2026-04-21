#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===============================================================================
PIPELINE METADATA
===============================================================================
*/
manifest {
    name        = 'wastewater-viral-pipeline'
    author      = 'BhardwajTulika'
    version     = '1.0.0'
    description = 'Metagenomic viral detection pipeline (HPC-ready)'
}

/*
===============================================================================
PARAMETERS (NO HARD-CODED HPC PATHS)
===============================================================================
*/
params.reads            = null
params.outdir           = "results"

params.kraken_db        = null
params.bbmap_refs       = null

params.genomad_db       = null
params.checkv_db        = null

params.megahit_min_len  = 500

/*
===============================================================================
VALIDATION
===============================================================================
*/
if (!params.reads) error "ERROR: --reads is required"
if (!params.kraken_db) error "ERROR: --kraken_db is required"
if (!params.bbmap_refs) error "ERROR: --bbmap_refs is required"

/*
===============================================================================
CHANNELS
===============================================================================
*/
read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
ref_files  = Channel.fromPath(params.bbmap_refs, checkIfExists: true)

/*
===============================================================================
FASTQC
===============================================================================
*/
process FASTQC {

    container "biocontainers/fastqc:v0.11.9_cv8"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("*_fastqc.zip"),
          path("*_fastqc.html")

    script:
    """
    set -euo pipefail
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}

/*
===============================================================================
FASTP
===============================================================================
*/
process FASTP {

    container "biocontainers/fastp:0.23.2"

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1_trimmed.fastq.gz"),
          path("${sample_id}_R2_trimmed.fastq.gz"),
          emit: trimmed_reads

    script:
    """
    set -euo pipefail

    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 ${sample_id}_R1_trimmed.fastq.gz \
        --out2 ${sample_id}_R2_trimmed.fastq.gz \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --html ${sample_id}.fastp.html \
        --json ${sample_id}.fastp.json
    """
}

/*
===============================================================================
RIBODETECTOR
===============================================================================
*/
process RIBODETECTOR {

    container "your_docker_or_singularity/ribodetector:0.3.1"

    publishDir "${params.outdir}/ribodetector", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.nonrrna.1.fq"),
          path("${sample_id}.nonrrna.2.fq"),
          emit: nonrrna_reads

    script:
    """
    ribodetector_cpu \
        -t ${task.cpus} \
        -l 100 \
        -i ${trimmed_reads[0]} ${trimmed_reads[1]} \
        -e both \
        --chunk_size 256 \
        -o ${sample_id}.nonrrna.1.fq ${sample_id}.nonrrna.2.fq \
        -r ${sample_id}.rrna.1.fq ${sample_id}.rrna.2.fq \
        --log ${sample_id}.log
    """
}

/*
===============================================================================
KRAKEN2
===============================================================================
*/
process KRAKEN2 {

    container "staphb/kraken2:2.1.3"

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}.report")
    path("${sample_id}.output")

    script:
    """
    set -euo pipefail

    kraken2 \
        --db ${params.kraken_db} \
        --paired ${reads[0]} ${reads[1]} \
        --threads ${task.cpus} \
        --report ${sample_id}.report \
        --output ${sample_id}.output
    """
}

/*
===============================================================================
BBMAP (TRUE PARALLEL PER REFERENCE)
===============================================================================
*/
process BBMAP {

    container "biocontainers/bbmap"

    publishDir "${params.outdir}/bbmap", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(ref_name), path(ref)

    output:
    path("${sample_id}.${ref_name}.covstats.txt")
    path("${sample_id}.${ref_name}.scafstats.txt")

    script:
    """
    set -euo pipefail

    bbmap.sh \
        ref=${ref} \
        in1=${reads[0]} \
        in2=${reads[1]} \
        threads=${task.cpus} \
        covstats=${sample_id}.${ref_name}.covstats.txt \
        scafstats=${sample_id}.${ref_name}.scafstats.txt
    """
}

/*
===============================================================================
MEGAHIT
===============================================================================
*/
process MEGAHIT {

    container "biocontainers/megahit:1.2.9"

    publishDir "${params.outdir}/megahit", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fa"), emit: contigs

    script:
    """
    set -euo pipefail

    megahit \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o ${sample_id}_asm \
        -t ${task.cpus} \
        --min-contig-len ${params.megahit_min_len}

    cp ${sample_id}_asm/final.contigs.fa ${sample_id}.contigs.fa
    """
}

/*
===============================================================================
GENOMAD
===============================================================================
*/
process GENOMAD {

    container "genomad:1.8.1"

    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), emit: viral

    script:
    """
    set -euo pipefail

    genomad end-to-end \
        ${contigs} \
        ${sample_id}_genomad \
        ${params.genomad_db} \
        --threads ${task.cpus}

    FILE=\$(find ${sample_id}_genomad -name "*virus.fna" | sort | head -n 1)

    cp \$FILE viral_contigs.fasta
    """
}

/*
===============================================================================
CHECKV
===============================================================================
*/
process CHECKV {

    container "biocontainers/checkv"

    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    tuple val(sample_id), path(viral)

    output:
    path("${sample_id}_checkv")

    script:
    """
    set -euo pipefail

    checkv end_to_end \
        ${viral} \
        ${sample_id}_checkv \
        -t ${task.cpus} \
        -d ${params.checkv_db}
    """
}

/*
===============================================================================
WORKFLOW (ARC / SLURM READY DAG)
===============================================================================
*/
workflow {

    read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ref_files  = Channel.fromPath(params.bbmap_refs)

    qc      = FASTQC(read_pairs)
    trimmed = FASTP(read_pairs)

    ribo    = RIBODETECTOR(trimmed.out.trimmed_reads)

    kraken  = KRAKEN2(ribo.out.nonrrna_reads)

    bbmap   = BBMAP(
                ribo.out.nonrrna_reads
                    .combine(ref_files)
                    .flatten()
              )

    asm     = MEGAHIT(ribo.out.nonrrna_reads)

    viral   = GENOMAD(asm.out.contigs)

    CHECKV(viral.out.viral)
}
