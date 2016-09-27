#!/usr/bin/env bash

time gzip -dc tracks/S1_LGL_MP.insert_length.gz \
  > tracks/S1_LGL_MP.insert_length.sparse.scratch

time genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = input tracks/S1_LGL_MP.insert_length.sparse.scratch --missing=-inf \
      = addconst 8000 \
      = slidingsum W=1000 D=W \
      = percentile 60 W=25 --min=1/inf --quiet \
      = input tracks/S1_LGL_MP.insert_length.sparse.scratch --missing=-inf \
      = addconst 8000 \
      = slidingsum W=1000 D=W \
      = percentile 3,5 W=25 --min=1/inf --quiet \
      = input tracks/S1_LGL_MP.insert_length.sparse.scratch --missing=-inf --destroy \
      = addconst 8000 \
      = clip --min=percentile3 \
      = mask tracks/hg19.Ns.dat --mask=percentile60 \
      = mask tracks/repeat_masker.hg19.dat --mask=percentile60 \
      = anticlump --average=percentile5 L=500 \
      = mask tracks/hg19.Ns.dat --mask=0 \
      = mask tracks/repeat_masker.hg19.dat --mask=0 \
  > tracks/S1_LGL_MP.insert_length.sparse.dat
