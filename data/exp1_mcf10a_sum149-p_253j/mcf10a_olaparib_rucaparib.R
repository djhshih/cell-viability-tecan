root <- "../.."

source(file.path(root, "R/common.R"));

cell_line <- "MCF10A";
compounds <- c("Olaparib", "Rucaparib");

data_t0_fname <- "tsv/mcf10a_d0_i1_2019-02-08.tsv";
data_t_fname <- "tsv/mcf10a_d5_i1_olaparib_rucaparib_2019-02-13.tsv";


date.format <- "%Y-%m-%d %H:%M";
time_t0 <- strptime("2019-02-08 16:00", format=date.format);
time_t <- strptime("2019-02-13 15:00", format=date.format);

assay <- list(
	method = "Presto Blue, fluorescence",
	incubation_hours = 1,
	duration_hours = as.numeric(difftime(time_t, time_t0, units="hours"))
);

annot_fname <- "annot/cell-viability-drug_2019-02-08b.csv";

design_t0_fname <- file.path(root, "design/plate-design_w66_blank-r.tsv");
design_t_fname <- file.path(root, "design/plate-design_w66_cp-tb_blank-r.tsv");

out_dir <- "norm";

ys <- normalize_viability(
	data_t0_fname, data_t_fname,
	design_t0_fname, design_t_fname,
	annot_fname,
	cell_line,
	compounds,
	assay
);

lapply(compounds, function(cp) {
		qwrite(ys[[cp]],
			filename("viability",
				path=out_dir,
				tag=c(tolower(cell_line), tolower(cp)),
				ext="rds"
			)
		);
	}
);

