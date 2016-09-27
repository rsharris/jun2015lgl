#!/usr/bin/env bash

time samtools view mapped/S1_LGL_MP.name_sorted.bam \
  | filtered_sam_to_intervals \
      --namesorted \
      --mergemates=300..6000 \
      --prohibit:"(CIGAR == *)" \
      --require:" (RNEXT == =)" \
      --require:" (PORIENT==T2T)" \
      --nonames \
      --progress=2M \
  | genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --novalue --show:uncovered \
  | tee tracks/S1_LGL_MP.short_inserts.depth.temp \
  | gzip \
  > tracks/S1_LGL_MP.short_inserts.depth.gz
