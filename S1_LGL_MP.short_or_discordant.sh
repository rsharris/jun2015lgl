#!/usr/bin/env bash

time cat tracks/S1_LGL_MP.short_inserts.depth.dense.dat \
  | genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = add tracks/S1_LGL_MP.discordant_mates.dense.dat \
      = binarize \
  > tracks/S1_LGL_MP.short_or_discordant.dat

