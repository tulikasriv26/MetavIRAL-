#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================
PARAMETERS (FROM CONFIG)
========================
*/

params.reads            = params.reads
params.outdir           = params.outdir
params.kraken_db        = params.kraken_db
params.bbmap_refs       = params.bbmap_refs
params.genomad_db       = params.genomad_db
params.checkv_db        = params.checkv_db
params.megahit_min_len  = params.megahit_min_len

/*
========================
CHANNELS
========================
*/

read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
ref_files  = Channel.fromPath(params.bbmap_refs, checkIfExists: true).collect()

/*
========================
FASTQC
========================
*/
process FASTQC {

    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    script:
    """
    set -euo pipefail
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}
/*
========================
FASTP
========================
*/
process FASTP {

    container "https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1_trimmed.fastq.gz"),
          path("${sample_id}_R2_trimmed.fastq.gz"),
          emit: trimmed_reads

    path "${sample_id}.fastp.html"
    path "${sample_id}.fastp.json"

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
} 
/*
========================
RIBODETECTOR
========================
*/
process RIBODETECTOR {

    container "https://depot.galaxyproject.org/singularity/ribodetector:0.3.1--pyhdfd78af_0"

    publishDir "${params.outdir}/ribodetector/non_rRNA", pattern: "*.nonrrna.*.fq", mode: 'copy'
    publishDir "${params.outdir}/ribodetector/rRNA",     pattern: "*.rrna.*.fq",    mode: 'copy'
    publishDir "${params.outdir}/ribodetector/logs",     pattern: "*.log",          mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}.nonrrna.1.fq"),
          path("${sample_id}.nonrrna.2.fq"),
          emit: nonrrna_reads

    tuple val(sample_id),
          path("${sample_id}.rrna.1.fq"),
          path("${sample_id}.rrna.2.fq"),
          emit: rrna_reads

    path "${sample_id}.log"

    script:
    """
    set -euo pipefail

    ribodetector_cpu \\
        -t ${task.cpus} \\
        -l 100 \\
        -i ${read1} ${read2} \\
        -e both \\
        --chunk_size 256 \\
        -o ${sample_id}.nonrrna.1.fq ${sample_id}.nonrrna.2.fq \\
        -r ${sample_id}.rrna.1.fq ${sample_id}.rrna.2.fq \\
        --log ${sample_id}.log
    """
}

/*
========================
KRAKEN2
========================
*/
process KRAKEN2_NONRRNA {

    container "https://depot.galaxyproject.org/singularity/kraken2:2.1.3--pl5321h077b44d_4"

    containerOptions "--bind /work:/work"

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${sample_id}.kraken2.report"
    path "${sample_id}.kraken2.output"

    script:
    """
    set -euo pipefail

    kraken2 \\
        --db ${params.kraken_db} \\
        --paired ${read1} ${read2} \\
        --threads ${task.cpus} \\
        --report ${sample_id}.kraken2.report \\
        --output ${sample_id}.kraken2.output
    """
}

/*
========================
BBMAP_NONRRNA
========================
*/
process BBMAP_NONRRNA {

    container "https://depot.galaxyproject.org/singularity/bbmap:39.13--he5f24ec_1"

    publishDir "${params.outdir}/bbmap", pattern: "*.sam.gz", mode: 'copy'
    publishDir "${params.outdir}/bbmap", pattern: "*.covstats.txt", mode: 'copy'
    publishDir "${params.outdir}/bbmap", pattern: "*.scafstats.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path refs

    output:
    path("${sample_id}.*.sam.gz")
    path "${sample_id}.*.covstats.txt"
    path "${sample_id}.*.scafstats.txt"

    script:
    def commands = refs.collect { ref ->
        def refName = ref.baseName
        """
        bbmap.sh \\
            ref=${ref} \\
            in1=${read1} \\
            in2=${read2} \\
            threads=${task.cpus} \\
            out=${sample_id}.${refName}.sam \\
            covstats=${sample_id}.${refName}.covstats.txt \\
            scafstats=${sample_id}.${refName}.scafstats.txt

        gzip -f ${sample_id}.${refName}.sam
        """
    }.join('\n')

    """
    set -euo pipefail
    ${commands}
    """
}
/*
========================
MEGAHIT
========================
*/
process MEGAHIT_NONRRNA {

    container "https://depot.galaxyproject.org/singularity/megahit:1.2.9--h43eeafb_5"

    publishDir "${params.outdir}/megahit", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("final.contigs.fa"),
          emit: contigs

    script:
    """
    set -euo pipefail

    megahit \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${sample_id}_assembly \\
        -t ${task.cpus} \\
        --presets meta-sensitive \\
        --min-contig-len ${params.megahit_min_len}

    cp ${sample_id}_assembly/final.contigs.fa final.contigs.fa
    """
}
/*
========================
GENOMAD
========================
*/
process GENOMAD_NONRRNA {

    container "https://depot.galaxyproject.org/singularity/genomad:1.8.1--pyhdfd78af_0"
    containerOptions "--bind /work:/work"
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), emit: viral_contigs

    script:
    """
    set -euo pipefail

    # Run geNomad
    genomad end-to-end \\
        ${contigs} \\
        ${sample_id}_genomad \\
        ${params.genomad_db} \\
        --threads ${task.cpus} \\
        --quiet

    # Extract viral contigs (more precise path)
    VIRAL_CONTIGS=\$(find ${sample_id}_genomad -path "*_summary/*_virus.fna" | head -n 1)

    if [ -z "\${VIRAL_CONTIGS}" ] || [ ! -f "\${VIRAL_CONTIGS}" ]; then
        echo "ERROR: No viral contigs found from geNomad" >&2
        exit 1
    fi

    cp "\${VIRAL_CONTIGS}" viral_contigs.fasta
    """
}
/*
========================
CHECKV
========================
*/
process CHECKV_NONRRNA {

    container "https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0"
    containerOptions "--bind /work:/work"
    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    tuple val(sample_id), path(viral_contigs)

    output:
    path("${sample_id}_checkv")

    script:
    """
    set -euo pipefail

    # Check DB exists (prevents your previous error)
    if [ ! -d "${params.checkv_db}" ]; then
        echo "ERROR: CheckV DB not found at ${params.checkv_db}" >&2
        exit 1
    fi

    checkv end_to_end \\
        ${viral_contigs} \\
        ${sample_id}_checkv \\
        -t ${task.cpus} \\
        -d ${params.checkv_db}
    """
}

/*
========================
WORKFLOW
========================
*/

workflow {

    FASTQC(read_pairs)
    FASTP(read_pairs)

    RIBODETECTOR(FASTP.out.trimmed_reads)

    KRAKEN2_NONRRNA(RIBODETECTOR.out.nonrrna_reads)
    BBMAP_NONRRNA(RIBODETECTOR.out.nonrrna_reads, ref_files)

    MEGAHIT_NONRRNA(RIBODETECTOR.out.nonrrna_reads)
    GENOMAD_NONRRNA(MEGAHIT_NONRRNA.out.contigs)
    CHECKV_NONRRNA(GENOMAD_NONRRNA.out.viral_contigs)
}
    """
}
