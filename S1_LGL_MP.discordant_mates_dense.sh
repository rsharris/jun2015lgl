#!/usr/bin/env bash

time genodsp \
      --chromosomes=data/hg19.chrom_lengths \
      --show:uncovered \
      = input discordant/S1_LGL_MP.BDB.MMQ40.MCP40.rmdup.bedgraph \
      = binarize \
      = slidingsum W=8000 D=W \
      = percentile 95 W=25 --min=1/inf --quiet \
      = input discordant/S1_LGL_MP.BDB.MMQ40.MCP40.rmdup.bedgraph \
      = binarize \
      = clump --average=percentile95 L=8000 \
  > tracks/S1_LGL_MP.discordant_mates.dense.dat
