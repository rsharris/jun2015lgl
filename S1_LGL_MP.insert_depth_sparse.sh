#!/usr/bin/env bash

time gzip -dc tracks/S1_LGL_MP.normal_inserts.depth.gz \
  > tracks/S1_LGL_MP.inserts.depth.sparse.scratch

time genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = input tracks/S1_LGL_MP.inserts.depth.sparse.scratch \
      = slidingsum W=1000 D=W \
      = percentile 60 W=25 --min=1/inf --quiet \
      = input tracks/S1_LGL_MP.inserts.depth.sparse.scratch \
      = slidingsum W=1000 D=W \
      = percentile 10 W=25 --min=1/inf --quiet \
      = input tracks/S1_LGL_MP.inserts.depth.sparse.scratch --destroy \
      = mask tracks/hg19.Ns.dat --mask=percentile60 \
      = mask tracks/repeat_masker.hg19.dat --mask=percentile60 \
      = anticlump --average=percentile10 L=500 \
      = mask tracks/hg19.Ns.dat --mask=0 \
      = mask tracks/repeat_masker.hg19.dat --mask=0 \
  > tracks/S1_LGL_MP.normal_inserts.depth.sparse.dat
