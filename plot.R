library(ggplot2)
library(drc)
library(io)
library(dplyr)
library(ggalt)
library(RColorBrewer)
library(digest)

hash_stamp <- function(x, trim=TRUE) {
	h <- sha1(x);
	if (trim) {
		strtrim(h, 7)
	} else {
		h
	}
}


ic50 <- function(fit) {
	ED(fit, 0.5, interval="delta", type="absolute")
}

ec50 <- function(fit) {
	ED(fit, 50, interval="delta", type="relative")
}

format_concentration <- function(x, unit = "uM") {
	if (x < 1) {
		x <- x * 1e3;
		# TODO make more robust
		unit <- "nM";
	}
	sprintf("%s %s", format(x, digits=3), unit)
}

fit_band <- function(fit, xlim) {
	d <- expand.grid(concentration = 10^seq(log10(xlim[1]), log10(xlim[2]), length=100));
	# TODO this generates many warnings
	pm <- suppressWarnings( predict(fit, newdata=d, interval="confidence") );
	data.frame(
		d,
		p = pm[,1],
		pmin = pm[,2],
		pmax = pm[,3]
	);
}


measure <- "relative_viability";
#measure <- "relative_growth";
#measure <- "relative_growth_rate";

draw.ci <- TRUE;

model <- "LL4";

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_talazoparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_selumetinib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i4/viability_sum149-p_selumetinib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i4/viability_sum149-p_selumetinib.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_selumetinib.rds"
);



#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_selumetinib.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_selumetinib.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_selumetinib.rds"
#);
#
#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-c2_selumetinib.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_selumetinib.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_selumetinib.rds"
#);
#
#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_vincristine.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vincristine.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_vincristine.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_vincristine.rds"
#);

#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_talazoparib.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_talazoparib.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_talazoparib.rds"
#);
#
#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_olaparib.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_olaparib.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_olaparib.rds"
#);
#
#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_cisplatin.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_cisplatin.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_cisplatin.rds"
#);

#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_vincristine.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_vincristine.rds"
#);

#fnames <- c(
#	"data/exp4_sum149-p_sum149-c2/norm/i4/viability_sum149-p_vincristine.rds",
#	"data/exp4_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_vincristine.rds"
#);

#fnames <- c(
#	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vincristine.rds",
#	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vincristine.rds"
#);


#fnames <- c(
#	i1 = "data/exp3_sum149-p_sum149-c2/norm/i1/viability_sum149-p_cisplatin.rds",
#	i2 = "data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
#	i3 = "data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_cisplatin.rds",
#	i4 = "data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_cisplatin.rds"
#);

#fnames <- c(
#	MCF10A_1="data/exp1_mcf10a_sum149-p_253j/norm/i1/viability_mcf10a_selumetinib.rds",
#	MCF10A_2="data/exp2_mcf10a_253j/norm/i1/viability_mcf10a_selumetinib.rds",
#	"SUM149-P"="data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_selumetinib.rds",
#	"SUM149-C2"="data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_selumetinib.rds"
#);
#
#fnames <- c(
#	exp1="data/exp1_mcf10a_sum149-p_253j/norm/i2/viability_mcf10a_cisplatin.rds",
#	exp2="data/exp2_mcf10a_253j/norm/i2/viability_mcf10a_cisplatin.rds"
#);
#
#fnames <- c(
#	exp1="data/exp1_mcf10a_sum149-p_253j/norm/i2/viability_mcf10a_vincristine.rds",
#	exp2="data/exp2_mcf10a_253j/norm/i2/viability_mcf10a_vincristine.rds"
#);

#fnames <- c(
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-p_cisplatin.rds",
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_cisplatin.rds"
#);
#
#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_cisplatin.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_cisplatin.rds"
#);
#
#fnames <- c(
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_rucaparib.rds",
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_rucaparib.rds"
#);
#
#fnames <- c(
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-p_rucaparib.rds",
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_rucaparib.rds"
#);
#
#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_rucaparib.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_rucaparib.rds"
#);
#
#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_olaparib.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_olaparib.rds"
#);
#
#fnames <- c(
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_olaparib.rds",
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_olaparib.rds"
#);
#
#fnames <- c(
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_talazoparib.rds",
#	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_talazoparib.rds"
#);
#
#fnames <- c(
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-p_talazoparib.rds",
#	"data/exp4_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_talazoparib.rds"
#);
#
#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_talazoparib.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_talazoparib.rds"
#);


#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_vincristine.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_talazoparib.rds"
#);
#
#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_vincristine.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_talazoparib.rds"
#);

#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_cisplatin.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_cisplatin.rds"
#);

#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_olaparib.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_olaparib.rds"
#);

#fnames <- c(
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-p_rucaparib.rds",
#	"data/exp3_sum149-p_sum149-c2/norm/i4/viability_sum149-c2_rucaparib.rds"
#);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_bms599626.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_bms599626.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_sch772984.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_sch772984.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vx11e.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vx11e.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_fr180204.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_fr180204.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_pictilisib.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_pictilisib.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_everolimus.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_everolimus.rds"
);


fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_cisplatin.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_cisplatin.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_vincristine.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_vincristine.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-p_talazoparib.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-p_olaparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i3/viability_sum149-c2_olaparib.rds"
);


fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_selumetinib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_talazoparib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_olaparib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_olaparib.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_olaparib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_olaparib.rds"
);


fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_cisplatin.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_cisplatin.rds"
);

fnames <- c(
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp5_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_cisplatin.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_cisplatin.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_everolimus.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_everolimus.rds"
);

fnames <- c(
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-p_everolimus.rds",
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_everolimus.rds"
);

fnames <- c(
	"data/exp11_sum149/norm/i2/viability_sum149-p_everolimus.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c2_everolimus.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_everolimus.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_everolimus.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_pictilisib.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_pictilisib.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-p_pictilisib.rds",
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_pictilisib.rds"
);

fnames <- c(
	"data/exp11_sum149/norm/i2/viability_sum149-p_pictilisib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c2_pictilisib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_pictilisib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_pictilisib.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_sch772984.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_sch772984.rds"
);

fnames <- c(
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-p_sch772984.rds",
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_sch772984.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp11_sum149/norm/i2/viability_sum149-p_sch772984.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c2_sch772984.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_sch772984.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_sch772984.rds"
);

fnames <- c(
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vx11e.rds",
	"data/exp8_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vx11e.rds"
);

fnames <- c(
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vx11e.rds",
	"data/exp9_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vx11e.rds"
);

fnames <- c(
	"data/exp11_sum149/norm/i2/viability_sum149-p_vx11e.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c2_vx11e.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_vx11e.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_vx11e.rds"
);

fnames <- c(
	"data/exp11_sum149/norm/i2/viability_sum149-p_lapatinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c2_lapatinib.rds"
);

fnames <- c(
	"data/exp12_sum149/norm/i2/viability_sum149-p_lapatinib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c2_lapatinib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_lapatinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_lapatinib.rds"
);


fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
	"data/exp10_sum149-c19_sum149-c30/norm/i2/viability_sum149-c19_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	"data/exp10_sum149-c19_sum149-c30/norm/i2/viability_sum149-c30_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_selumetinib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c19_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_selumetinib.rds"
);

fnames <- c(
	# FIXME exp5 cannot be integered due to different number of columns!
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c19_selumetinib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c30_selumetinib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_selumetinib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c19_talazoparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_talazoparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_talazoparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_talazoparib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c19_talazoparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c20_talazoparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c29_talazoparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c30_talazoparib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_talazoparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_talazoparib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_olaparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c19_olaparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_olaparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_olaparib.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_olaparib.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_olaparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c19_olaparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c20_olaparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c29_olaparib.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c30_olaparib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_olaparib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_olaparib.rds"
);

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_cisplatin.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c19_cisplatin.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_cisplatin.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_cisplatin.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_cisplatin.rds"
);

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_cisplatin.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c19_cisplatin.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c20_cisplatin.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c29_cisplatin.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c30_cisplatin.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_cisplatin.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_cisplatin.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vincristine.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vincristine.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c19_vincristine.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c20_vincristine.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c29_vincristine.rds",
	"data/exp11_sum149/norm/i2/viability_sum149-c30_vincristine.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vincristine.rds",
	"data/exp4_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vincristine.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c19_vincristine.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c20_vincristine.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c29_vincristine.rds",
	"data/exp12_sum149/norm/i2/viability_sum149-c30_vincristine.rds"
);
#draw.ci <- FALSE;

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_vincristine.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_vincristine.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_vincristine.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_vincristine.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_vincristine.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_vincristine.rds"
);
#draw.ci <- FALSE;


fnames <- c(
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds",
	#"data/exp11_sum149/norm/i2/viability_sum149-c19_talazoparib.rds",
	#"data/exp11_sum149/norm/i2/viability_sum149-c20_talazoparib.rds",
	#"data/exp11_sum149/norm/i2/viability_sum149-c29_talazoparib.rds",
	#"data/exp11_sum149/norm/i2/viability_sum149-c30_talazoparib.rds"
	"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q2-pre-t_talazoparib.rds",
	#"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q4-pre-t_talazoparib.rds",
	"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q6-pre-t_talazoparib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_selumetinib.rds",
	#"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	#"data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
	#"data/exp13_sum149/norm/i2/viability_sum149-c19_selumetinib.rds",
	#"data/exp13_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	#"data/exp13_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	#"data/exp13_sum149/norm/i2/viability_sum149-c30_selumetinib.rds",
	"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q2-pre-t_selumetinib.rds",
	#"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q4-pre-t_selumetinib.rds",
	"data/exp15_sum149-q-pre/norm/i2/viability_sum149-q6-pre-t_selumetinib.rds"
);

fnames <- c(
	"data/exp13_sum149/norm/i2/viability_sum149-p_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c2_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c19_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c20_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c29_selumetinib.rds",
	"data/exp13_sum149/norm/i2/viability_sum149-c30_selumetinib.rds"
);



xs <- lapply(fnames, qread);

samples <- unlist(lapply(xs, function(x) x$meta$cell_line));
compounds <- unlist(lapply(xs, function(x) x$meta$compound));

if (length(unique(compounds)) == 1) {
	pl.title <- compounds[1];
	if (is.null(names(xs))) {
		names(xs) <- samples;
	}
	group.name <- "sample";
} else if (length(unique(samples)) == 1) {
	pl.title <- samples[1];
	if (is.null(names(xs))) {
		names(xs) <- compounds;
	}
	group.name <- "compound";
} else {
	pl.title <- NULL;
}

sources <- lapply(xs, function(x) x$meta$sources);

out.fname <- filename("drc",
	tag = c(
		gsub("_", "-", measure),
		paste0(tolower(samples), "_", tolower(compounds)),
		hash_stamp(sources)
	)
);

# select measure
xs <- lapply(xs, function(x) { x$data$y <- x$data[[measure]]; x });

# fix the upper asymptote (d = 1)
if (model == "LL5") {
	fixed <- c(NA, NA, 1, NA, NA);
	rv.fits <- lapply(xs,
		function(x) {
			drm(y ~ concentration, data = x$data, fct = LL.5(fixed=fixed))
		}
	);
} else {
	# fix the upper asymptote (d = 1)
	fixed <- c(NA, NA, 1, NA);
	rv.fits <- lapply(xs,
		function(x) {
			drm(y ~ concentration, data = x$data, fct = LL.4(fixed=fixed))
		}
	);
}

ic50s <- data.frame(group = names(xs), do.call(rbind, lapply(rv.fits, ic50)), row.names=NULL);
ec50s <- data.frame(group = names(xs), do.call(rbind, lapply(rv.fits, ec50)), row.names=NULL);

concentrations <- unlist(lapply(xs, function(x) x$data$concentration));
max.conc <- max(concentrations);
min.conc <- min(concentrations);

bands <- lapply(rv.fits, function(fit) fit_band(fit, c(min.conc, max.conc)));

options(scipen=1e4)

stopifnot(length(xs) <= 7);
cols <- brewer.pal(7, "Set1")[1:length(xs)];
names(cols) <- names(xs);


all.n <- do.call(rbind,
	mapply(function(x, g) data.frame(group = g, x$data), xs, names(xs), SIMPLIFY=FALSE)
);
rownames(all.n) <- NULL;

z <- qnorm(0.975);
all.n.s <- group_by(all.n, group, concentration) %>%
	summarize(ym = mean(y), yse = sd(y) / length(y), ymin = ym + z*yse, ymax = ym - z*yse);


ymin <- min(all.n.s$ym, all.n.s$ymin, unlist(lapply(bands, function(b) c(b$p, b$pmin))));
ymax <- max(all.n.s$ym, all.n.s$ymax, unlist(lapply(bands, function(b) c(b$p, b$pmax))));


g <- ggplot(all.n.s, aes(x = concentration, y = ym)) + theme_bw() +
		geom_point(aes(colour=group)) +
		scale_x_log10() + 
		scale_fill_manual(values=cols) +
		scale_colour_manual(values=cols) + 
		guides(fill = FALSE) +
		ylab(gsub("_", " ", tools::toTitleCase(measure))) +
		xlab("Concentration (uM)") +
		ggtitle(pl.title);

if (draw.ci) {
	data.ymin <- min(ymin);
	data.ymax <- max(ymax);
} else {
	data.ymin <- min(all.n.s$ym);
	data.ymax <- max(all.n.s$ym);
}
	
if (measure == "relative_viability") {
	g <- g + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
		limits=c(min(0, data.ymin), max(1.0, data.ymax)));
} else {
	g <- g + scale_y_continuous(breaks=c(-1.0, -0.5, 0, 0.5, 1.0),
		limits=c(min(-1.0, data.ymin), max(1.0, data.ymax)));
}

for (i in 1:length(xs)) {
	g <- g + geom_line(data=bands[[i]], aes(x = concentration, y = p), colour = cols[i], size=0.75);
 	if (draw.ci) {
		g <- g + geom_ribbon(data=bands[[i]], aes(x = concentration, y = p, ymin = pmin, ymax = pmax), fill = cols[i], alpha=0.2);
	}
	ic50 <- ic50s$Estimate[i];
	if (!is.na(ic50) && ic50 <= max.conc && ic50 >= min.conc) {
		g <- g + geom_vline(xintercept=ic50, colour=cols[i], alpha=0.5) +
			annotate("text", x = ic50 * 0.9, y = 0.05 * (i - 1), hjust=1, label=format_concentration(ic50));
	}
}

qdraw(g, width = 6, height = 5, file = insert(out.fname, ext="pdf"));

qwrite(ic50s, insert(out.fname, "ic50", ext="tsv"));
qwrite(ec50s, insert(out.fname, "ec50", ext="tsv"));


ggplot(all.n, aes(x=group, y=value)) +
	geom_point() + theme_bw()

