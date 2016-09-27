#!/usr/bin/env bash

time cat tracks/S1_LGL_MP.insert_length_sparse_or_normal_inserts_sparse.dat \
  | proximal_feature_intervals --positive \
      tracks/S1_LGL_MP.short_or_discordant.dat \
      --proximity=2K \
  | awk '{ print $1,$2,$3 }' \
  | keep_first \
  | close_intervals 2K \
  | fill_genomic_interval_gaps --chroms=data/hg19.chrom_lengths \
  > tracks/S1.called_insertions.dat

specs="IL-:tracks/S1_LGL_MP.insert_length.sparse.dat
       NI-:tracks/S1_LGL_MP.normal_inserts.depth.sparse.dat
       SI+:tracks/S1_LGL_MP.short_inserts.depth.dense.dat
       SM+:tracks/S1_LGL_MP.discordant_mates.dense.dat"
echo ${specs} \
  | tr " :" "\n " \
  | while read nick trackname ; do
      cat tracks/S1.called_insertions.dat \
        | proximal_feature_intervals --positive \
            ${trackname} \
            --proximity=2K \
        | awk '{ print $1"~"$2"~"$3,nick,1 }' nick=${nick} \
        | keep_first
      done \
  | collect_tags --separator=~ \
  | awk '{ print $1,"#",$2 }' \
  | sed "s/~/ /g" \
  | encodachrom | env LC_ALL=C sort -k 1,1n -k 2,2n | decodachrom \
  | intervals_to_ucsc_catalog \
      --show:numbers \
      --show:comments \
      --genome=hg19 \
      --split=20,80 \
      --center=20K \
      --center=150K \
      --catalogonly \
      --title="called insertions (${today})" \
      tracks/S1.called_insertions.catalog.html
