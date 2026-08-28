options(stringsAsFactors = FALSE, warn = 1)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else normalizePath(".")
validation_dir <- normalizePath(file.path(script_dir), mustWork = TRUE)
output_dir <- Sys.getenv(
  "WICHERS_FIGURE_DIR",
  file.path(dirname(validation_dir), "outputs", "manuscript_figures")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(validation_dir, "data", "raw", "wichers_transition", "ESMdata.zip")
model_path <- file.path(validation_dir, "results", "wichers_transition_models.csv")
summary_path <- file.path(validation_dir, "results", "wichers_transition_summary.csv")

d <- read.csv(unz(zip_path, "ESMdata/ESMdata.csv"))
d$measurement_date <- as.Date(d$date, "%d/%m/%y")
d$experiment_day <- as.integer(d$measurement_date - min(d$measurement_date)) + 1
daily <- aggregate(
  rep(1L, nrow(d)),
  by = list(measurement_date = d$measurement_date, experiment_day = d$experiment_day),
  FUN = sum
)
names(daily)[3] <- "momentary_assessments"
daily$concentrat <- vapply(
  daily$measurement_date,
  function(day) {
    values <- d$concentrat[d$measurement_date == day]
    finite <- values[!is.na(values)]
    if (length(finite)) finite[1] else NA_real_
  },
  numeric(1)
)
daily <- daily[order(daily$experiment_day), ]
weekly <- unique(d[!is.na(d$dep), c("measurement_date", "experiment_day", "concentrat", "dep")])
weekly <- weekly[order(weekly$experiment_day), ]
models <- read.csv(model_path)
transition_summary <- read.csv(summary_path)

accent <- "#176B87"
ink <- "#27343A"
grey <- "#8A9499"
light <- "#DCECEF"
boundary <- transition_summary$best_boundary_day[1]
pre_mean <- transition_summary$pre_boundary_mean[1]
post_mean <- transition_summary$post_boundary_mean[1]
bootstrap_share <- transition_summary$bootstrap_proportion_day_127[1]

draw_figure <- function(path, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(path, width = 2100, height = 2250, res = 300, pointsize = 9, bg = "white")
  } else {
    pdf(path, width = 7.0, height = 7.5, pointsize = 9, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  layout(matrix(1:3, nrow = 3), heights = c(1.04, 1.18, 0.98))
  par(
    family = "sans",
    fg = ink,
    col.axis = ink,
    col.lab = ink,
    mar = c(3.8, 4.8, 2.8, 4.8),
    mgp = c(2.35, 0.65, 0),
    tcl = -0.25,
    cex.axis = 0.92,
    cex.lab = 1.0,
    font.axis = 1,
    font.lab = 1,
    lend = "round"
  )

  panel_heading <- function(label, title) {
    mtext(
      paste0("(", label, ")  ", title),
      side = 3,
      adj = 0,
      line = 0.65,
      font = 2,
      cex = 1.08
    )
  }

  barplot(
    daily$momentary_assessments,
    names.arg = rep("", nrow(daily)),
    border = NA,
    col = light,
    space = 0,
    axes = FALSE,
    xlab = "Study day",
    ylab = "Momentary assessments"
  )
  axis(2, at = c(0, 5, 10))
  day_to_x <- function(day) (day - min(daily$experiment_day)) / diff(range(daily$experiment_day)) * nrow(daily)
  axis(1, at = day_to_x(c(1, 60, 127, 180, 239)), labels = c(1, 60, 127, 180, 239))
  abline(v = day_to_x(boundary), lty = 2, lwd = 1.6, col = accent)
  par(new = TRUE)
  plot(
    day_to_x(daily$experiment_day),
    daily$concentrat,
    type = "s",
    lwd = 2.1,
    col = grey,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = c(0, 160)
  )
  axis(4, at = c(0, 50, 100, 150), col.axis = grey)
  mtext("Venlafaxine dose (mg)", side = 4, line = 2.75, col = grey, cex = 0.84)
  panel_heading("a", "Intensive observation and medication taper")
  legend("topright", legend = c("Assessments/day", "Medication dose", "Day 127"), fill = c(light, NA, NA), border = c(NA, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2.1, 1.6), col = c(light, grey, accent), bty = "n", cex = 0.76)

  par(mar = c(3.8, 4.8, 2.8, 4.8))
  plot(
    weekly$experiment_day,
    weekly$dep,
    type = "o",
    pch = 21,
    bg = "white",
    col = ink,
    lwd = 1.6,
    xlim = c(1, 239),
    ylim = range(c(weekly$dep, pre_mean, post_mean)) + c(-0.12, 0.12),
    xlab = "Study day",
    ylab = "Weekly depression score"
  )
  abline(v = boundary, lty = 2, lwd = 1.8, col = accent)
  segments(min(weekly$experiment_day), pre_mean, boundary, pre_mean, col = grey, lwd = 2.5)
  segments(boundary, post_mean, max(weekly$experiment_day), post_mean, col = accent, lwd = 2.5)
  text(
    boundary + 4,
    max(weekly$dep),
    sprintf("Published and recovered weekly occasion\nday 127; %.1f%% of bootstraps", 100 * bootstrap_share),
    adj = c(0, 1),
    col = accent,
    cex = 0.76
  )
  panel_heading("b", "The search selects the same weekly occasion")

  labels <- c(
    Step127 = "Day-127 step (search-adjusted)",
    AR_Step127 = "AR + day-127 step (adjusted)",
    DelayedDose = "Delayed discontinuation",
    AR_DelayedDose = "AR + delayed discontinuation",
    LinearTime = "Linear time",
    AR_Dose = "AR + medication dose",
    AR = "Autoregressive",
    LinearDose = "Linear medication dose",
    Constant = "Constant"
  )
  models$display <- unname(labels[models$model])
  models <- models[order(models$bic_adjusted, decreasing = TRUE), ]
  par(mar = c(4.0, 11.2, 2.8, 1.2))
  mids <- barplot(
    models$bic_adjusted,
    names.arg = models$display,
    horiz = TRUE,
    las = 1,
    xlab = "BIC (lower is better)",
    col = ifelse(models$model == "Step127", accent, grey),
    border = NA,
    xlim = c(0, max(models$bic_adjusted) * 1.13),
    cex.names = 0.77
  )
  text(
    models$bic_adjusted + 0.8,
    mids,
    sprintf("%.1f", models$bic_adjusted),
    adj = 0,
    cex = 0.74,
    col = ink
  )
  panel_heading("c", "The step is modestly favoured after the search penalty")
}

draw_figure(file.path(output_dir, "Fig_6_Wichers.png"), "png")
draw_figure(file.path(output_dir, "Fig_6_Wichers.pdf"), "pdf")
cat("Wichers manuscript figure written to", output_dir, "\n")
