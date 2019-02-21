#!/bin/bash

indir=$1
outdir=$indir/..

if (( $# < 1 )); then
	echo "usage: $0 <indir>" >&2
	exit -1
fi

mkdir -p $outdir/{csv,tsv}

for f in $indir/*.xlsx; do
	fname=${f##*/}
	fstem=${fname%.*}
	csv=$outdir/csv/${fstem}.csv
	if [[ ! -f $csv ]]; then
		echo $f
		xlsx2csv $f $csv
	fi
	tsv=$outdir/tsv/${fstem}.tsv
	if [[ ! -f $tsv ]]; then
		echo $csv
		R/tecan-convert.r $csv --outdir $outdir/tsv
	fi
done

