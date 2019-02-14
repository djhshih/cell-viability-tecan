
#' Compute apparent growth rate and doubling time
#' based on the basic exponential model.
# @param x   population at time t
# @param x0  population at time 0
# @param dt  elapsed time
growth_characteristics <- function(x, x0, dt) {
	r <- log(x / x0) / dt;
	list(growth_rate = r, doubling_time = log(2) / r)
}

#' Compute relative growth rate value.
# @param x   measurement at time t
# @param x0  measurement at time 0
# @param xc  control measurement at time t
gr_value <- function(x, x0, xc) {
	2^( log2(x / x0) / log2(xc / x0) ) - 1
}

#' Compute relative growth value.
#' @param x   measurement at time t
#' @param x0  measurement at time 0
#' @param xc  control measurement at time t
gi_value <- function(x, x0, xc) {
	ifelse( x > x0,
		(x - x0) / (xc - x0),
		(x - x0) / x0
	)
}

#' Normalize viability assay.
normalize_viability <- function(
	data_t0_fname, data_t_fname,
	design_t0_fname, design_t_fname,
	annot_fname,
	cell_line,
	compounds,
	assay
) {

annot <- qread(annot_fname);

annotf <- select(annot[match(compounds, annot$drug), ],
	compound = drug, max_conc, dilution_factor, solvent) %>%
	mutate(compound = factor(compound));

design_t0 <- qread(design_t0_fname);
design_t <- qread(design_t_fname);

d0 <- qread(data_t0_fname) %>% left_join(design_t0) %>% filter(!is.na(group));
dx <- qread(data_t_fname) %>% left_join(design_t) %>% filter(!is.na(group));

print(design_t)
print(dx)

d0.blank <- mean(filter(d0, group == "ctl_blank")$value);

d0 <- mutate(d0, value_bg = pmax(value - d0.blank, 0));
d0.untrt <- mean(filter(d0, group == "ctl_untrt")$value_bg);


dx$compound <- factor(dx$compound, labels=compounds);

dx.blank <- mean(filter(dx, group == "ctl_blank")$value);

dx <- mutate(dx, value_bg = pmax(value - dx.blank, 0));

dx.vehicle <- mean(filter(dx, group == "ctl_vehicle")$value_bg);
dx.untrt <- mean(filter(dx, group == "ctl_untrt")$value_bg);

growth <- growth_characteristics(dx.untrt, d0.untrt, assay$duration_hours);

dx <- mutate(dx,
	relative_viability = value_bg / dx.vehicle,
	relative_growth = gi_value(value_bg, d0.untrt, dx.vehicle),
	relative_growth_rate = gr_value(value_bg, d0.untrt, dx.vehicle)
);

dx.cps <- split(dx, dx$compound);
annotf.cps <- split(annotf, annotf$compound);

dilution_to_conc <- function(dilution, max_conc, dilution_factor) {
	max_conc * (1 / dilution_factor)^(dilution - 1)
}

ys <- lapply(compounds,
	function(cp) {

		dx.cps[[cp]]$concentration <- dilution_to_conc(
			dx.cps[[cp]]$dilution, annotf.cps[[cp]]$max_conc, annotf.cps[[cp]]$dilution_factor);

		list(
			meta = list(
				cell_line = cell_line,
				compound = cp,
				assay = assay,
				sources = list(
					design_t0_fname = design_t0_fname,
					design_t_fname = design_t_fname,
					data_t0_fname = data_t0_fname,
					data_t_fname = data_t_fname
				),
				growth = growth,
				ctl_values = list(
					t0 = list(
						blank = d0.blank,
						untrt = d0.untrt
					),
					t = list(
						blank = dx.blank,
						untrt = dx.untrt,
						vehicle = dx.vehicle
					)
				)
			),
			data = select(dx.cps[[cp]], -group, -dilution, -compound)
		);
	}
);

ys
}

