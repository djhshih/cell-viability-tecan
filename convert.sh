#!/bin/bash

outdir=tsv

mkdir -p $outdir

for f in csv/*.csv; do
	R/tecan-convert.r $f --outdir $outdir
done
