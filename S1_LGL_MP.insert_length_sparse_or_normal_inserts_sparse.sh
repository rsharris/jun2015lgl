#!/usr/bin/env bash

time cat tracks/S1_LGL_MP.insert_length.sparse.dat \
  | genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = add tracks/S1_LGL_MP.normal_inserts.depth.sparse.dat \
      = binarize \
  > tracks/S1_LGL_MP.insert_length_sparse_or_normal_inserts_sparse.dat
