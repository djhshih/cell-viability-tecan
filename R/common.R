library(io);
library(dplyr);

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

standardize_compound <- function(x) {
	ifelse(
		grepl("\\d", x),
		# name contains a digit; probably a preliminary compound id:
		#   use uppercase
		toupper(x),
		# name only contains letters: it is probably a generic compound name:
		#   use title case
		tools::toTitleCase(x)
	)
}

dilution_to_conc <- function(dilution, max_conc, dilution_factor) {
	max_conc * (1 / dilution_factor)^(dilution - 1)
}

#' Normalize viability assay.
normalize_viability <- function(
	data_t0_fname, data_t_fname,
	design_t0_fname, design_t_fname,
	annot_fname,
	cell_line,
	compounds,
	assay,
	cells = NULL
) {

annot <- qread(annot_fname, stringsAsFactors=TRUE);

annotf <- select(annot[match(compounds, annot$drug), ],
	compound = drug, max_conc, dilution_factor, solvent) %>%
	mutate(compound = factor(compound));

if (!is.null(design_t0_fname)) {
	design_t0 <- qread(design_t0_fname, stringsAsFactors=TRUE);
}
design_t <- qread(design_t_fname, stringsAsFactors=TRUE);

if (!is.null(data_t0_fname)) {
	d0 <- qread(data_t0_fname) %>% left_join(design_t0) %>% filter(!is.na(group));
} else {
	d0 <- NULL;
}
dx <- qread(data_t_fname) %>% left_join(design_t) %>% filter(!is.na(group));

if (!is.null(cells)) {
	if (!is.null(d0) && !is.null(d0$cell)) {
		d0 <- mutate(d0, cell = factor(cell, levels=levels(cell), labels=cells)) %>%
			filter(is.na(cell) | cell == cell_line) %>% select(-cell);
	}
	if (!is.null(dx$cell)) {
		dx <- mutate(dx, cell = factor(cell, levels=levels(cell), labels=cells)) %>%
			filter(is.na(cell) | cell == cell_line) %>% select(-cell);
	}
}

print(dx)

if (!is.null(d0)) {
	d0.blank <- mean(filter(d0, group == "ctl_blank")$value);

	d0 <- mutate(d0, value_bg = pmax(value - d0.blank, 0));
	d0.untrt <- mean(filter(d0, group == "ctl_untrt")$value_bg);
} else {
	d0.blank <- NULL;
	d0.untrt <- NULL;
}

# specifcy levels s.t. function does not throw error when dx$compound does not
# contain all the possible levels.
dx$compound <- factor(dx$compound, levels=levels(dx$compound), labels=compounds);

dx.blank <- mean(filter(dx, group == "ctl_blank")$value);

dx <- mutate(dx, value_bg = pmax(value - dx.blank, 0));

dx.vehicle <- mean(filter(dx, group == "ctl_vehicle")$value_bg);
dx.untrt <- mean(filter(dx, group == "ctl_untrt")$value_bg);

if (!is.null(d0.untrt)) {
	growth <- growth_characteristics(dx.untrt, d0.untrt, assay$duration_hours);
} else {
	growth <- NULL;
}

dx <- mutate(dx,
	relative_viability = value_bg / dx.vehicle,
);

if (!is.null(d0.untrt)) {
	dx <- mutate(dx,
		relative_growth = gi_value(value_bg, d0.untrt, dx.vehicle),
		relative_growth_rate = gr_value(value_bg, d0.untrt, dx.vehicle)
	);
}

dx.cps <- split(dx, dx$compound);
annotf.cps <- split(annotf, annotf$compound);

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
names(ys) <- compounds;

ys
}

#' Normalize viability assay for drug combinations.
normalize_viability_combo <- function(
	data_t0_fname, data_t_fname,
	design_t0_fname, design_t_fname,
	annot_fname,
	cell_line,
	combos,
	assay,
	cells = NULL
) {

annot <- qread(annot_fname, stringsAsFactors=TRUE);

combo_to_info <- function(combos) {
	ss <- strsplit(combos, "_", fixed=TRUE);
	combo_methods <- unlist(lapply(ss, function(x) x[1]));
	combo_rests <- unlist(lapply(ss, function(x) x[2]));
	combo_infos <- mapply(
		function(method, combo) {
			compounds <- strsplit(combo, "-", fixed=TRUE)[[1]];
			list(method = method, compounds = compounds)
		},
		combo_methods,
		combo_rests,
		SIMPLIFY = FALSE
	);
	names(combo_infos) <- combos;

	combo_infos
}

annotf <- select(annot, compound = drug, max_conc, dilution_factor, solvent);

if (!is.null(design_t0_fname)) {
	design_t0 <- qread(design_t0_fname, stringsAsFactors=TRUE);
}
design_t <- qread(design_t_fname, stringsAsFactors=TRUE);

if (!is.null(data_t0_fname)) {
	d0 <- qread(data_t0_fname) %>% left_join(design_t0) %>% filter(!is.na(group));
} else {
	d0 <- NULL;
}

dx <- qread(data_t_fname) %>% left_join(design_t) %>% filter(!is.na(group));

if (!is.null(cells)) {
	if (!is.null(d0) && !is.null(d0$cell)) {
		d0 <- mutate(d0, cell = factor(cell, levels=levels(cell), labels=cells)) %>%
			filter(is.na(cell) | cell == cell_line) %>% select(-cell);
	}
	if (!is.null(dx$cell)) {
		dx <- mutate(dx, cell = factor(cell, levels=levels(cell), labels=cells)) %>%
			filter(is.na(cell) | cell == cell_line) %>% select(-cell);
	}
}

print(dx)

if (!is.null(d0)) {
	d0.blank <- mean(filter(d0, group == "ctl_blank")$value);

	d0 <- mutate(d0, value_bg = pmax(value - d0.blank, 0));
	d0.untrt <- mean(filter(d0, group == "ctl_untrt")$value_bg);
} else {
	d0.blank <- NULL;
	d0.untrt <- NULL;
}

# specifcy levels s.t. function does not throw error when dx$compound does not
# contain all the possible levels.
# consider each combo as a new compound
dx$compound <- factor(dx$compound, levels=levels(dx$compound), labels=combos);

dx.blank <- mean(filter(dx, group == "ctl_blank")$value);

dx <- mutate(dx, value_bg = pmax(value - dx.blank, 0));

dx.vehicle <- mean(filter(dx, group == "ctl_vehicle")$value_bg);
dx.untrt <- mean(filter(dx, group == "ctl_untrt")$value_bg);

if (!is.null(d0.untrt)) {
	growth <- growth_characteristics(dx.untrt, d0.untrt, assay$duration_hours);
} else {
	growth <- NULL;
}

dx <- mutate(dx,
	relative_viability = value_bg / dx.vehicle,
);

if (!is.null(d0.untrt)) {
	dx <- mutate(dx,
		relative_growth = gi_value(value_bg, d0.untrt, dx.vehicle),
		relative_growth_rate = gr_value(value_bg, d0.untrt, dx.vehicle)
	);
}

dx.cps <- split(dx, dx$compound);
annotf.cps <- split(annotf, annotf$compound);


ys <- lapply(combos,
	function(combo) {

		info <- combo_to_info(combo)[[1]];

		for (i in 1:length(info$compounds)) {
			cp <- info$compounds[i];
			dx.cps[[combo]][[paste0("concentration", i)]] <- dilution_to_conc(
				dx.cps[[combo]][[paste0("dilution", i)]], annotf.cps[[cp]]$max_conc, annotf.cps[[cp]]$dilution_factor);
		}

		list(
			meta = list(
				cell_line = cell_line,
				compound = combo,
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
			data = select(dx.cps[[combo]], -group, -dilution1, -dilution2, -compound)
		);
	}
);
names(ys) <- combos;

ys
}
