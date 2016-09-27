#!/usr/bin/env bash

time gzip -dc tracks/S1_LGL_MP.short_inserts.depth.gz \
  > tracks/S1_LGL_MP.inserts.depth.dense.scratch

time genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = input tracks/S1_LGL_MP.inserts.depth.dense.scratch \
      = slidingsum W=1000 D=W \
      = percentile 90,92 W=25 --min=1/inf --quiet \
      = input tracks/S1_LGL_MP.inserts.depth.dense.scratch --destroy \
      = clip --max=percentile92 \
      = clump --average=percentile90 L=1000 \
  > tracks/S1_LGL_MP.short_inserts.depth.dense.dat
