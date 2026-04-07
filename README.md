# MetavIRAL-
This repository contains a Nextflow pipeline for preprocessing paired-end metatranscriptomics reads and running downstream analysis on the non-rRNA fraction.
Pipeline Overview

This workflow performs:

Quality control → FastQC
Read trimming → fastp
rRNA removal → RiboDetector
Taxonomic classification → Kraken2
Reference mapping → BBMap
Assembly → MEGAHIT
Viral identification → geNomad
Viral genome quality → CheckV
