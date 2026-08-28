options(stringsAsFactors = FALSE, warn = 1)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else normalizePath(".")
write_outputs <- Sys.getenv("VALIDATION_WRITE", "1") != "0"
result_dir <- Sys.getenv("VALIDATION_OUTPUT_DIR", file.path(root, "results"))

clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
entropy <- function(p) {
  p <- clip(p)
  -p * log(p) - (1 - p) * log(1 - p)
}
logloss <- function(y, p) -mean(y * log(clip(p)) + (1 - y) * log(clip(1 - p)))
brier <- function(y, p) mean((y - p)^2)
roc_auc <- function(y, p) {
  r <- rank(p, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
pr_auc <- function(y, p) {
  o <- order(p, decreasing = TRUE)
  yy <- y[o]
  if (!sum(yy)) return(NA_real_)
  tp <- cumsum(yy)
  recall <- tp / sum(yy)
  precision <- tp / seq_along(tp)
  sum((recall - c(0, head(recall, -1))) * precision)
}
cluster_delta <- function(x, p_new, p_base, seed = 240826L, reps = 10000L) {
  ids <- unique(x$id)
  d <- sapply(ids, function(who) {
    q <- x$id == who
    logloss(x$y[q], p_new[q]) - logloss(x$y[q], p_base[q])
  })
  set.seed(seed)
  boot <- replicate(reps, mean(sample(d, replace = TRUE)))
  c(
    mean = mean(d),
    low = unname(quantile(boot, .025)),
    high = unname(quantile(boot, .975)),
    better = sum(d < 0),
    participants = length(d)
  )
}

fit_ai <- function(tr, type) {
  if (type == "AI_exact") {
    fn <- function(par) {
      p <- plogis(-exp(par[1]) * (tr$risk_basis * par[2] + tr$ambiguity))
      -sum(tr$y * log(clip(p)) + (1 - tr$y) * log(clip(1 - p)))
    }
    fit <- optim(c(0, 0), fn, method = "BFGS", control = list(maxit = 1000))
  } else if (type == "AI_field") {
    fn <- function(par) {
      p <- plogis(par[1] - exp(par[2]) * (tr$risk_basis * par[3] + tr$ambiguity))
      -sum(tr$y * log(clip(p)) + (1 - tr$y) * log(clip(1 - p)))
    }
    fit <- optim(c(qlogis(mean(tr$y)), 0, 0), fn, method = "BFGS", control = list(maxit = 1000))
  } else {
    fn <- function(par) {
      p <- plogis(par[1] - exp(par[2]) * (tr$risk_basis * par[3] + tr$ambiguity) + par[4] * tr$prev)
      -sum(tr$y * log(clip(p)) + (1 - tr$y) * log(clip(1 - p)))
    }
    fit <- optim(c(qlogis(mean(tr$y)), 0, 0, 0), fn, method = "BFGS", control = list(maxit = 1000))
  }
  list(par = fit$par, type = type)
}
predict_ai <- function(fit, te) {
  p <- fit$par
  gamma <- exp(p[if (fit$type == "AI_exact") 1 else 2])
  dc <- p[if (fit$type == "AI_exact") 2 else 3]
  eta <- -gamma * (te$risk_basis * dc + te$ambiguity)
  if (fit$type != "AI_exact") eta <- p[1] + eta
  if (fit$type == "AI_field_persistence") eta <- eta + p[4] * te$prev
  plogis(eta)
}
lopo_action <- function(x, dataset, prior_total) {
  models <- c(
    "Intercept", "Markov", "ActionHistory", "OutcomeLearning",
    "AI_exact", "AI_field", "AI_field_persistence"
  )
  pred <- setNames(lapply(models, function(z) rep(NA_real_, nrow(x))), models)
  for (who in unique(x$id)) {
    tr <- x[x$id != who, ]
    te <- x[x$id == who, ]
    ii <- which(x$id == who)
    fits <- list(
      Intercept = glm(y ~ 1, tr, family = binomial()),
      Markov = glm(y ~ prev + log_gap, tr, family = binomial()),
      ActionHistory = glm(y ~ prev + action_order + log_gap, tr, family = binomial()),
      OutcomeLearning = glm(y ~ prev + good_difference + ambiguity + log_gap, tr, family = binomial())
    )
    for (m in names(fits)) pred[[m]][ii] <- predict(fits[[m]], te, type = "response")
    for (m in c("AI_exact", "AI_field", "AI_field_persistence")) {
      pred[[m]][ii] <- predict_ai(fit_ai(tr, m), te)
    }
  }
  stable <- x$y == x$prev
  delta <- cluster_delta(x, pred$AI_field_persistence, pred$Markov)
  delta_history <- cluster_delta(
    x,
    pred$AI_field_persistence,
    pred$ActionHistory,
    seed = 240827L
  )
  data.frame(
    dataset = dataset,
    prior_total = prior_total,
    model = models,
    n = nrow(x),
    participants = length(unique(x$id)),
    action_prevalence = mean(x$y),
    changes = sum(!stable),
    log_loss = sapply(models, function(m) logloss(x$y, pred[[m]])),
    brier = sapply(models, function(m) brier(x$y, pred[[m]])),
    change_log_loss = sapply(models, function(m) logloss(x$y[!stable], pred[[m]][!stable])),
    stable_log_loss = sapply(models, function(m) logloss(x$y[stable], pred[[m]][stable])),
    ai_minus_markov_cluster_mean = delta["mean"],
    ai_minus_markov_ci_low = delta["low"],
    ai_minus_markov_ci_high = delta["high"],
    participants_ai_better = delta["better"],
    ai_minus_action_history_cluster_mean = delta_history["mean"],
    ai_minus_action_history_ci_low = delta_history["low"],
    ai_minus_action_history_ci_high = delta_history["high"],
    participants_ai_better_than_action_history = delta_history["better"]
  )
}

make_fisher <- function(d, prior_total = 2, outcome_name = "content", cutoff = 50) {
  d$timestamp <- as.POSIXct(d$start, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  d$action <- as.integer(rowMeans(d[c("avoid_activity", "avoid_people")], na.rm = FALSE) < cutoff)
  d$good <- as.integer(d[[outcome_name]] >= cutoff)
  d <- d[complete.cases(d[c("id", "timestamp", "action", "good")]), ]
  half <- prior_total / 2
  out <- list()
  k <- 1
  for (who in unique(d$id)) {
    g <- d[d$id == who, ]
    g <- g[order(g$timestamp), ]
    eg <- eb <- wb <- wg <- half
    ce <- cw <- 0
    for (i in seq_len(nrow(g))) {
      a <- eg / (eg + eb)
      b <- wb / (wb + wg)
      if (i > 1) {
        z <- if (ce + cw) (ce - cw) / (ce + cw) else 0
        gap <- as.numeric(difftime(g$timestamp[i], g$timestamp[i - 1], units = "hours"))
        out[[k]] <- data.frame(
          id = who,
          y = g$action[i],
          prev = g$action[i - 1],
          action_order = z,
          a = a,
          b = b,
          risk_basis = 1 - a - b,
          ambiguity = entropy(a) - entropy(b),
          good_difference = a - (1 - b),
          log_gap = log1p(pmax(gap, 0))
        )
        k <- k + 1
      }
      if (g$action[i] == 1) {
        if (g$good[i] == 1) eg <- eg + 1 else eb <- eb + 1
        ce <- ce + 1
      } else {
        if (g$good[i] == 0) wb <- wb + 1 else wg <- wg + 1
        cw <- cw + 1
      }
    }
  }
  do.call(rbind, out)
}

make_jang_action <- function(d, prior_total = 2) {
  half <- prior_total / 2
  out <- list()
  k <- 1
  for (who in unique(d$id)) {
    g <- d[d$id == who, ]
    g <- g[order(g$date), ]
    eg <- eb <- wb <- wg <- half
    ce <- cw <- 0
    for (i in seq_len(nrow(g))) {
      a <- eg / (eg + eb)
      b <- wb / (wb + wg)
      consecutive <- i > 1 && as.integer(g$date[i] - g$date[i - 1]) == 1
      if (consecutive) {
        z <- if (ce + cw) (ce - cw) / (ce + cw) else 0
        out[[k]] <- data.frame(
          id = who,
          y = g$exercise[i],
          prev = g$exercise[i - 1],
          action_order = z,
          a = a,
          b = b,
          risk_basis = 1 - a - b,
          ambiguity = entropy(a) - entropy(b),
          good_difference = a - (1 - b),
          log_gap = log(2)
        )
        k <- k + 1
      }
      if (consecutive) {
        previous <- g$exercise[i - 1]
        safe <- 1 - g$panic[i]
        if (previous == 1) {
          if (safe == 1) eg <- eg + 1 else eb <- eb + 1
        } else {
          if (safe == 0) wb <- wb + 1 else wg <- wg + 1
        }
      }
      if (g$exercise[i] == 1) ce <- ce + 1 else cw <- cw + 1
    }
  }
  do.call(rbind, out)
}

make_jang_event <- function(d, prior = 1) {
  out <- list()
  k <- 1
  for (who in unique(d$id)) {
    g <- d[d$id == who, ]
    g <- g[order(g$date), ]
    eg <- eb <- wb <- wg <- panic_yes <- panic_no <- prior
    ce <- cw <- 0
    for (i in seq_len(nrow(g))) {
      consecutive <- i > 1 && as.integer(g$date[i] - g$date[i - 1]) == 1
      if (consecutive) {
        previous <- g$exercise[i - 1]
        safe <- 1 - g$panic[i]
        if (previous == 1) {
          if (safe == 1) eg <- eg + 1 else eb <- eb + 1
        } else {
          if (safe == 0) wb <- wb + 1 else wg <- wg + 1
        }
      }
      panic_yes <- panic_yes + g$panic[i]
      panic_no <- panic_no + 1 - g$panic[i]
      a <- eg / (eg + eb)
      b <- wb / (wb + wg)
      z <- if (ce + cw) (ce - cw) / (ce + cw) else 0
      next_day <- i < nrow(g) && as.integer(g$date[i + 1] - g$date[i]) == 1
      if (next_day) {
        out[[k]] <- data.frame(
          id = who,
          y = g$panic[i + 1],
          panic = g$panic[i],
          exercise = g$exercise[i],
          anxiety = g$anxiety[i],
          negative = g$negative_feeling[i],
          positive = g$positive_feeling[i],
          log_steps = log1p(g$steps_mean[i]),
          action_order = z,
          mechanism_risk = if (g$exercise[i] == 1) 1 - a else b,
          online_risk = panic_yes / (panic_yes + panic_no)
        )
        k <- k + 1
      }
      if (g$exercise[i] == 1) ce <- ce + 1 else cw <- cw + 1
    }
  }
  do.call(rbind, out)
}

lopo_jang_event <- function(x) {
  x <- x[complete.cases(x), ]
  x$lm <- qlogis(clip(x$mechanism_risk))
  x$lo <- qlogis(clip(x$online_risk))
  models <- c(
    "Intercept", "PanicMarkov", "Symptoms", "Activity", "ActionHistory",
    "MechanismDirect", "OnlineDirect", "CalibratedMechanism",
    "CalibratedOnline", "OnlinePlusMechanism", "Full"
  )
  pred <- setNames(lapply(models, function(z) rep(NA_real_, nrow(x))), models)
  for (who in unique(x$id)) {
    tr <- x[x$id != who, ]
    te <- x[x$id == who, ]
    ii <- which(x$id == who)
    fits <- list(
      Intercept = glm(y ~ 1, tr, family = binomial()),
      PanicMarkov = glm(y ~ panic, tr, family = binomial()),
      Symptoms = glm(y ~ panic + anxiety + negative + positive, tr, family = binomial()),
      Activity = glm(y ~ panic + exercise + log_steps, tr, family = binomial()),
      ActionHistory = glm(y ~ panic + action_order, tr, family = binomial()),
      CalibratedMechanism = glm(y ~ panic + lm, tr, family = binomial()),
      CalibratedOnline = glm(y ~ panic + lo, tr, family = binomial()),
      OnlinePlusMechanism = glm(y ~ panic + lo + lm, tr, family = binomial()),
      Full = glm(
        y ~ panic + anxiety + negative + positive + exercise + log_steps + action_order + lm,
        tr,
        family = binomial()
      )
    )
    for (m in names(fits)) pred[[m]][ii] <- predict(fits[[m]], te, type = "response")
    pred$MechanismDirect[ii] <- te$mechanism_risk
    pred$OnlineDirect[ii] <- te$online_risk
  }
  d1 <- cluster_delta(x, pred$CalibratedMechanism, pred$CalibratedOnline)
  d2 <- cluster_delta(x, pred$OnlinePlusMechanism, pred$CalibratedOnline)
  data.frame(
    dataset = "Jang_2024_next_day_panic",
    model = models,
    n = nrow(x),
    participants = length(unique(x$id)),
    events = sum(x$y),
    event_prevalence = mean(x$y),
    log_loss = sapply(models, function(m) logloss(x$y, pred[[m]])),
    brier = sapply(models, function(m) brier(x$y, pred[[m]])),
    roc_auc = sapply(models, function(m) roc_auc(x$y, pred[[m]])),
    pr_auc = sapply(models, function(m) pr_auc(x$y, pred[[m]])),
    mechanism_minus_online_cluster_mean = d1["mean"],
    mechanism_minus_online_ci_low = d1["low"],
    mechanism_minus_online_ci_high = d1["high"],
    online_plus_mechanism_minus_online_cluster_mean = d2["mean"],
    online_plus_mechanism_minus_online_ci_low = d2["low"],
    online_plus_mechanism_minus_online_ci_high = d2["high"]
  )
}

wichers_transition <- function(zip_path) {
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
      value <- d$concentrat[d$measurement_date == day]
      finite <- value[!is.na(value)]
      if (length(finite)) finite[1] else NA_real_
    },
    numeric(1)
  )
  daily <- daily[order(daily$experiment_day), ]
  w <- unique(d[!is.na(d$dep), c("measurement_date", "experiment_day", "concentrat", "dep")])
  w <- w[order(w$experiment_day), ]
  w$post_127 <- as.integer(w$experiment_day > 127)
  w$lag <- c(NA, head(w$dep, -1))
  w$gap <- c(NA, diff(w$experiment_day))
  zero_dose_day <- min(w$experiment_day[w$concentrat == 0], na.rm = TRUE)
  w$log_days_since_zero_dose <- log1p(pmax(w$experiment_day - zero_dose_day, 0))
  common <- w[complete.cases(w), ]
  fits <- list(
    Constant = lm(dep ~ 1, common),
    LinearDose = lm(dep ~ I(concentrat / 150), common),
    LinearTime = lm(dep ~ experiment_day, common),
    DelayedDose = lm(dep ~ log_days_since_zero_dose, common),
    Step127 = lm(dep ~ post_127, common),
    AR = lm(dep ~ lag + gap, common),
    AR_Dose = lm(dep ~ lag + gap + I(concentrat / 150), common),
    AR_DelayedDose = lm(dep ~ lag + gap + log_days_since_zero_dose, common),
    AR_Step127 = lm(dep ~ lag + gap + post_127, common)
  )
  searched_breakpoint <- grepl("Step127", names(fits))
  aicc_raw <- sapply(fits, function(f) {
    n_fit <- nobs(f)
    parameters <- length(coef(f)) + 1
    AIC(f) + 2 * parameters * (parameters + 1) / (n_fit - parameters - 1)
  })
  bic_raw <- sapply(fits, BIC)
  models <- data.frame(
    dataset = "Wichers_2016_N_of_1",
    model = names(fits),
    n = sapply(fits, nobs),
    rmse = sapply(fits, function(f) sqrt(mean(residuals(f)^2))),
    aicc_raw = aicc_raw,
    aicc_adjusted = aicc_raw + ifelse(searched_breakpoint, 2, 0),
    bic_raw = bic_raw,
    bic_adjusted = bic_raw + ifelse(searched_breakpoint, log(nrow(common)), 0),
    searched_breakpoint_parameters = as.integer(searched_breakpoint),
    r_squared = sapply(fits, function(f) summary(f)$r.squared)
  )
  n <- nrow(w)
  candidates <- 4:(n - 4)
  split_sse <- function(values, cut) {
    sum((values[1:cut] - mean(values[1:cut]))^2) +
      sum((values[(cut + 1):length(values)] - mean(values[(cut + 1):length(values)]))^2)
  }
  max_f <- function(values) {
    residual_null <- sum((values - mean(values))^2)
    residual_split <- sapply(candidates, function(cut) split_sse(values, cut))
    ((residual_null - residual_split) / 1) / (residual_split / (length(values) - 2))
  }
  sse <- sapply(candidates, function(cut) split_sse(w$dep, cut))
  k <- candidates[which.min(sse)]
  raw_f_profile <- max_f(w$dep)
  profile <- data.frame(
    boundary_day = w$experiment_day[candidates],
    residual_sum_of_squares = sse,
    max_f_statistic = raw_f_profile
  )
  set.seed(1272016L)
  permutation_reps <- 10000L
  raw_max_f <- max(raw_f_profile)
  raw_permuted <- replicate(permutation_reps, max(max_f(sample(w$dep))))
  raw_p <- (1 + sum(raw_permuted >= raw_max_f)) / (permutation_reps + 1)

  detrended <- residuals(lm(dep ~ experiment_day, w))
  detrended_max_f <- max(max_f(detrended))
  detrended_permuted <- replicate(
    permutation_reps,
    max(max_f(sample(detrended)))
  )
  detrended_p <- (1 + sum(detrended_permuted >= detrended_max_f)) /
    (permutation_reps + 1)

  fitted_step <- c(
    rep(mean(w$dep[1:k]), k),
    rep(mean(w$dep[(k + 1):n]), n - k)
  )
  step_residuals <- w$dep - fitted_step
  step_residuals <- step_residuals - mean(step_residuals)
  bootstrap_reps <- 10000L
  bootstrap_days <- replicate(bootstrap_reps, {
    simulated <- fitted_step + sample(step_residuals, replace = TRUE)
    selected_cut <- candidates[which.min(sapply(
      candidates,
      function(cut) split_sse(simulated, cut)
    ))]
    w$experiment_day[selected_cut]
  })
  bootstrap <- as.data.frame(table(bootstrap_days), stringsAsFactors = FALSE)
  names(bootstrap) <- c("boundary_day", "bootstrap_count")
  bootstrap$boundary_day <- as.integer(bootstrap$boundary_day)
  bootstrap$bootstrap_proportion <- bootstrap$bootstrap_count / bootstrap_reps
  bootstrap <- bootstrap[order(-bootstrap$bootstrap_count, bootstrap$boundary_day), ]
  coverage_index <- which(cumsum(bootstrap$bootstrap_proportion) >= .95)[1]
  highest_frequency_95_set <- sort(bootstrap$boundary_day[seq_len(coverage_index)])
  bootstrap_day_127 <- bootstrap$bootstrap_proportion[bootstrap$boundary_day == 127]

  pre <- w$dep[w$concentrat == 0 & w$experiment_day <= 127]
  post <- w$dep[w$concentrat == 0 & w$experiment_day > 127]
  summary <- data.frame(
    momentary_assessments = nrow(d),
    observed_days = nrow(daily),
    best_boundary_day = w$experiment_day[k],
    best_boundary_date = as.character(w$measurement_date[k]),
    selected_weekly_occasion = k,
    paper_boundary_day = 127,
    weekly_observations = nrow(w),
    common_model_observations = nrow(common),
    zero_dose_day = zero_dose_day,
    pre_boundary_mean = mean(w$dep[1:k]),
    post_boundary_mean = mean(w$dep[(k + 1):n]),
    same_zero_dose_pre_n = length(pre),
    same_zero_dose_post_n = length(post),
    same_zero_dose_pre_mean = mean(pre),
    same_zero_dose_post_mean = mean(post),
    same_zero_dose_difference = mean(post) - mean(pre),
    step127_bic_raw_common_n = bic_raw["Step127"],
    step127_bic_adjusted_common_n = models$bic_adjusted[models$model == "Step127"],
    delayed_dose_bic_common_n = bic_raw["DelayedDose"],
    delta_bic_delayed_minus_adjusted_step = bic_raw["DelayedDose"] -
      models$bic_adjusted[models$model == "Step127"],
    raw_max_f = raw_max_f,
    raw_max_f_permutation_p = raw_p,
    detrended_max_f = detrended_max_f,
    detrended_max_f_permutation_p = detrended_p,
    bootstrap_reps = bootstrap_reps,
    bootstrap_proportion_day_127 = bootstrap_day_127,
    bootstrap_highest_frequency_95_set = paste(highest_frequency_95_set, collapse = ";")
  )
  inference <- data.frame(
    test = c("max_F_against_constant_mean", "max_F_after_linear_detrending"),
    statistic = c(raw_max_f, detrended_max_f),
    permutation_reps = permutation_reps,
    p_value = c(raw_p, detrended_p)
  )
  list(
    models = models,
    summary = summary,
    inference = inference,
    bootstrap = bootstrap,
    daily = daily,
    weekly = w,
    profile = profile
  )
}

fisher_raw <- read.delim(
  file.path(root, "data/raw/fisher_openesm/0033_fisher_ts.tsv"),
  check.names = FALSE
)
fisher_results <- do.call(rbind, lapply(c(2, 10, 40), function(prior) {
  lopo_action(
    make_fisher(fisher_raw, prior),
    "Fisher_2017_avoidance_contentment",
    prior
  )
}))

jang_raw <- read.delim(
  file.path(root, "data/raw/jang_openesm/0017_jang_ts.tsv"),
  check.names = FALSE
)
jang_raw <- jang_raw[complete.cases(jang_raw[c("id", "date", "exercise", "panic")]), ]
jang_raw$date <- as.Date(jang_raw$date)
jang_action_results <- do.call(rbind, lapply(c(2, 10, 40), function(prior) {
  lopo_action(
    make_jang_action(jang_raw, prior),
    "Jang_2024_exercise_next_day_safety",
    prior
  )
}))

jang_event_results <- lopo_jang_event(make_jang_event(jang_raw))

wichers <- wichers_transition(
  file.path(root, "data/raw/wichers_transition/ESMdata.zip")
)

if (write_outputs) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(fisher_results, file.path(result_dir, "fisher_action_models.csv"), row.names = FALSE)
  write.csv(jang_action_results, file.path(result_dir, "jang_action_models.csv"), row.names = FALSE)
  write.csv(jang_event_results, file.path(result_dir, "jang_event_models.csv"), row.names = FALSE)
  write.csv(wichers$models, file.path(result_dir, "wichers_transition_models.csv"), row.names = FALSE)
  write.csv(wichers$summary, file.path(result_dir, "wichers_transition_summary.csv"), row.names = FALSE)
  write.csv(wichers$inference, file.path(result_dir, "wichers_change_point_inference.csv"), row.names = FALSE)
  write.csv(wichers$bootstrap, file.path(result_dir, "wichers_boundary_bootstrap.csv"), row.names = FALSE)
  write.csv(wichers$daily, file.path(result_dir, "wichers_daily_density.csv"), row.names = FALSE)
  write.csv(wichers$weekly, file.path(result_dir, "wichers_weekly_series.csv"), row.names = FALSE)
  write.csv(wichers$profile, file.path(result_dir, "wichers_change_point_profile.csv"), row.names = FALSE)
}

cat("\nFisher action models\n")
print(fisher_results[order(fisher_results$prior_total, fisher_results$log_loss), ], row.names = FALSE)
cat("\nJang action models\n")
print(jang_action_results[order(jang_action_results$prior_total, jang_action_results$log_loss), ], row.names = FALSE)
cat("\nJang next-day panic models\n")
print(jang_event_results[order(jang_event_results$log_loss), ], row.names = FALSE)
cat("\nWichers transition models\n")
print(wichers$models[order(wichers$models$bic_adjusted), ], row.names = FALSE)
print(wichers$summary, row.names = FALSE)
