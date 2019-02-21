#!/usr/bin/env Rscript

library(argparser);
library(filenamer);

pr <- arg_parser("Convert Tecan matrix data csv file to well list tsv file.");
pr <- add_argument(pr, "input", help = "input csv file");
pr <- add_argument(pr, "--outdir", help = "output directory", default=".");

argv <- parse_args(pr);

in.fname <- argv$input;

out.fname <- as.filename(in.fname);
out.fname$ext <- "tsv";
out.fname$path <- argv$outdir;

dir.create(out.fname$path, showWarnings=FALSE, recursive=TRUE);

x <- readLines(in.fname);

# extract matrix block
first.idx <- grep("^<>", x);
last.idx <- grep("(^,)|(^$)", x[first.idx:length(x)])[1] + first.idx - 2;
xf <- x[first.idx:last.idx];

sep <- ",";
tokens <- strsplit(xf, sep);

cols <- tokens[[1]][-1];
rows <- unlist(lapply(tokens[-1], function(x) x[1]));
d <- lapply(tokens[-1], function(x) x[-1]);

wells <- paste0(rep(rows, each=length(cols)), rep(cols, times=length(rows)));

y <- data.frame(
	well = wells,
	value = as.numeric(unlist(d))
);

write.table(y, tag(out.fname), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE);

