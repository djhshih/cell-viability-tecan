#!/bin/bash

indir=$1
outdir=$indir/../tsv

if (( $# < 1 )); then
	echo "usage: $0 <indir>" >&2
	exit -1
fi

mkdir -p $outdir

for f in $indir/*.csv; do
	R/tecan-convert.r $f --outdir $outdir
done
