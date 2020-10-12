library(ggplot2)
library(drc)
library(io)
library(dplyr)
library(ggalt)
library(RColorBrewer)


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

model <- "LL4";

# third fname must point to the results from the combination

fnames <- c(
	"Talazoparib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"Selumetinib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_selumetinib.rds",
	"Talazoparib-Selumetinib"="data/exp6_sum149-p_sum149-c2/norm/i2/viability_sum149-p_combo-con_talazoparib-selumetinib.rds"
);

fnames <- c(
	"Talazoparib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds",
	"Selumetinib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_selumetinib.rds",
	"Talazoparib-Selumetinib"="data/exp6_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_combo-con_talazoparib-selumetinib.rds"
);

fnames <- c(
	"Talazoparib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_talazoparib.rds",
	"Vincristine"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-p_vincristine.rds",
	"Talazoparib-Vincristine"="data/exp6_sum149-p_sum149-c2/norm/i2/viability_sum149-p_combo-con_talazoparib-vincristine.rds"
);

fnames <- c(
	"Talazoparib"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_talazoparib.rds",
	"Vincristine"="data/exp3_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_vincristine.rds",
	"Talazoparib-Vincristine"="data/exp6_sum149-p_sum149-c2/norm/i2/viability_sum149-c2_combo-con_talazoparib-vincristine.rds"
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

out.fname <- filename("drc", tag=c(gsub("_", "-", measure), paste0(tolower(samples), "_", tolower(compounds))));

x.combo <- xs[[3]];
concentrations1 <- unique(sort(x.combo$data$concentration1));
concentrations2 <- unique(sort(x.combo$data$concentration2));

ratio <- concentrations1 / concentrations2;

# we only support the case where the two drugs are used a constant ratio
stopifnot(abs(ratio - ratio[1]) < 1e-6)
ratio <- ratio[1];

# convert concentration2 to concentration1
xs[[2]]$data$concentration <- xs[[2]]$data$concentration * ratio;
# use concentration1
xs[[3]]$data$concentration <- xs[[3]]$data$concentration1;
xs[[3]]$data$concentration1 <- NULL;
xs[[3]]$data$concentration2 <- NULL;

# select measure
xs <- lapply(xs, function(x) {
	x$data$y <- x$data[[measure]];
	x
});

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

concentrations <- concentrations1;
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
		scale_fill_manual(values=cols) +
		scale_colour_manual(values=cols) + 
		guides(fill = FALSE) +
		ylab(gsub("_", " ", tools::toTitleCase(measure))) +
		ggtitle(pl.title) + 
		theme(legend.position = "bottom", legend.title = element_blank()) +
		scale_x_log10(
			name = sprintf("%s concentration (uM)", names(xs)[1]),
			sec.axis = sec_axis(
				trans =  ~ . / ratio,
				name = sprintf("%s concentration (uM)", names(xs)[2])
			)
		)
	
if (measure == "relative_viability") {
	g <- g + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
		limits=c(min(0, ymin), max(1.0, ymax)));
} else {
	g <- g + scale_y_continuous(breaks=c(-1.0, -0.5, 0, 0.5, 1.0),
		limits=c(min(-1.0, ymin), max(1.0, ymax)));
}

for (i in 1:length(xs)) {
	g <- g + geom_line(data=bands[[i]], aes(x = concentration, y = p), colour = cols[i]) +
		geom_ribbon(data=bands[[i]], aes(x = concentration, y = p, ymin = pmin, ymax = pmax), fill = cols[i], alpha=0.2);
	ic50 <- ic50s$Estimate[i];
	if (!is.na(ic50) && ic50 <= max.conc && ic50 >= min.conc) {
		g <- g + geom_vline(xintercept=ic50, colour=cols[i], alpha=0.5) +
			annotate("text", x = ic50 * 0.9, y = 0.05 * (i - 1), hjust=1, label=format_concentration(ic50));
	}
}

qdraw(g, width = 5, height = 6, file = insert(out.fname, ext="pdf"));

qwrite(ic50s, insert(out.fname, "ic50", ext="tsv"));
qwrite(ec50s, insert(out.fname, "ec50", ext="tsv"));


ggplot(all.n, aes(x=group, y=value)) +
	geom_point() + theme_bw()

