library(io)

root <- "../../";

source(file.path(root, "R/common.R"));


# global parameters

cells <- c("MCF10A", "253J", "SUM149-P");
date.format <- "%Y-%m-%d %H:%M";
time_t0 <- strptime("2019-02-08 16:00", format=date.format);
time_t <- strptime("2019-02-13 15:00", format=date.format);
indir <- "tsv";
annot_fname <- "annot/cell-viability-drug_2019-02-08b.csv";
design_t0_fname <- file.path(root, "design/plate-design_w66_blank-r.tsv");
design_t_fname <- file.path(root, "design/plate-design_w66_cp-tb_blank-r.tsv");
assay <- list(
	method = "Presto Blue, fluorescence",
	incubation_hours = NA,
	duration_hours = as.numeric(difftime(time_t, time_t0, units="hours"))
);
pattern_t <- "d5";
pattern_t0 <- "d0";

# enumerate input files

in.fnames <- list_files(indir);

in.fnames.t <- grep(pattern_t, in.fnames, value=TRUE);
in.fnames.t0 <- grep(pattern_t0, in.fnames, value=TRUE);

parse_fname <- function(x) {
	lapply(
		strsplit(x, "_", fixed=TRUE),
		function(tokens) {
			p <- list(
				cell_line = toupper(tokens[1]),
				day = tokens[2],
				incubation = tokens[3],
				compounds = standardize_compound(tokens[4:(length(tokens)-1)]),
				assay = assay,
				cells = cells
			);
			p$assay$incubation_hours = as.integer(sub("i", "", tokens[3], fixed=TRUE));

			p
		}
	)
}

params <- parse_fname(in.fnames.t);

params <- mapply(
	function(fname, p) {
		p$data_t_fname <- fname;
		p$data_t0_fname <- grep(p$incubation, in.fnames.t0, value=TRUE)[1];
		p
	},
	in.fnames.t,
	params,
	SIMPLIFY=FALSE
);


for (p in params) {
	print(p)

	out_dir <- file.path("norm", p$incubation);

	# normalize each compound separately
	ys <- normalize_viability(
		file.path(indir, p$data_t0_fname), file.path(indir, p$data_t_fname),
		design_t0_fname, design_t_fname,
		annot_fname,
		p$cell_line,
		p$compounds,
		p$assay,
		cells = p$cells
	);

	# write each compound separately
	lapply(p$compounds, function(cp) {
			qwrite(ys[[cp]],
				filename("viability",
					path = out_dir,
					tag = c(tolower(p$cell_line), tolower(cp)),
					ext = "rds",
					date = NA
				)
			);
		}
	);

}

