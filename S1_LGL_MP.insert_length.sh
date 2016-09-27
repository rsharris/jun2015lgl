#!/usr/bin/env bash

chrSizes=`cat data/hg19.chrom_lengths | awk '{ print $1":"$2 }'`

time samtools view /afs/bx.psu.edu/depot/data/poss_lab/LGL_new/mapped/S1_LGL_MP.name_sorted.bam \
  | filtered_sam_to_intervals \
      --namesorted \
      --mergemates=300..30000 \
      --prohibit:"(CIGAR == *)" \
      --require:" (RNEXT == =)" \
      --require:" (PORIENT==T2T)" \
      --nonames \
      --progress=2M \
  | awk '{ print $0,$3-$2 }' \
  | chrom_avg L=0 ${chrSizes} --precision=0 \
  | awk '{ print $1,$2,$3,($4-8000) }' \
  | tr " " "\t" \
  | tee tracks/S1_LGL_MP.insert_length.temp \
  | gzip \
  > tracks/S1_LGL_MP.insert_length.gz
