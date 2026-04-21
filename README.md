# MetavIRAL-
This repository contains a Nextflow pipeline for preprocessing paired-end metatranscriptomics reads and running downstream analysis on the non-rRNA fraction.

# Pipeline Overview
This workflow performs
1.Quality control → FastQC
2.Read trimming → fastp
3.rRNA, non-rRNA identification → RiboDetector
4.Taxonomic classification → Kraken2 
5.Reference mapping → BBMap
6.Assembly → MEGAHIT
7.Viral identification → geNomad
8.Viral genome quality → CheckV

# Input Requirements:
1. Paired-end FASTQ files:
*_R1_001.fastq.gz
*_R2_001.fastq.gz

3. Reference genomes (for BBMap):
*.fasta (e.g., 8 viral genomes)

5. Databases:
(a) Kraken2 database (https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz)
(b) geNomad database (https://zenodo.org/records/14886553)
(c) CheckV database (https://portal.nersc.gov/CheckV/)

# Parameters
| Parameter                | Description                                 |
| ------------------------ | ------------------------------------------- |
| `params.reads`           | Input FASTQ/FASTQ.GZ pattern                |
| `params.outdir`          | Output directory                            |
| `params.sif_dir`         | Directory containing Singularity containers |
| `params.kraken_db`       | Kraken2 database path                       |
| `params.bbmap_refs`      | Reference genomes (FASTA)                   |
| `params.genomad_db`      | geNomad database                            |
| `params.checkv_db`       | CheckV database                             |
| `params.megahit_min_len` | Minimum contig length                       |

### Required containers:

- fastqc: biocontainers/fastqc:v0.11.9_cv8  
- fastp: biocontainers/fastp:0.23.2  
- kraken2: staphb/kraken2:2.1.3  
- megahit: biocontainers/megahit:1.2.9  
- bbmap: biocontainers/bbmap 
-All tools are executed via Singularity container



