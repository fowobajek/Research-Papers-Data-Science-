# SAP3: Mediation of the cohort postponement of first births by girls' schooling
# Interventional effects (Vansteelandt & Daniel 2017) on restricted mean age at first birth (to age 30)
source(file.path(Sys.getenv("SAP_ROOT", "."), "R", "00_config.R"))
suppressPackageStartupMessages({ library(rstpm2); library(survival); library(nnet); library(parallel) })
set.seed(SEED); RNGkind("L'Ecuyer-CMRG")
o <- file.path(OUT, "sap3"); wcsv <- function(x, f) write_csv(x, file.path(o, f))
if (B_SAP3 < 200L) stop("BOOT_B_SAP3 must be at least 200; do not report percentile confidence intervals from a small bootstrap.")
ir <- readRDS(file.path(OUT, "data", "sap3_women.rds")) |>
  filter(!is.na(cohort), !is.na(sec), !is.na(zone), !is.na(religion), !is.na(t_fb)) |>
  mutate(zone = factor(zone, ZONES), religion = factor(religion),
         # time since age 10, administratively censored at age 30 (the RMST horizon)
         event = as.integer(event == 1 & t_fb <= 30), time = pmin(t_fb, 30) - 10)
wcsv(count(ir, cohort, year, name = "women"), "sap3_analytic_sample.csv")

# ---- Descriptive: cumulative probability of first birth by 15, 18, 20 by cohort (weighted KM) ----
km <- survfit(Surv(t_fb, event) ~ cohort, data = ir, weights = wt)
s <- summary(km, times = c(15, 18, 20), extend = TRUE)
cum <- tibble(cohort = sub("cohort=", "", s$strata), age = s$time, pct = 100 * (1 - s$surv),
              lo = 100 * (1 - s$upper), hi = 100 * (1 - s$lower), at_risk = s$n.risk) |>
  mutate(pct = ifelse(at_risk < 25, NA, pct))   # suppress where too few women observed to that age
wcsv(cum, "sap3_table1_cumulative_first_birth.csv")

# ---- Mediation machinery ----
TMAX <- 20; grid <- seq(0, TMAX, by = .25)   # time since age 10 -> RMST to age 30
d2 <- filter(ir, cohort %in% c("1950-1959", "1990-1999")) |> mutate(A = as.integer(cohort == "1990-1999"))
rmst_profiles <- function(fit, nd) {
  long <- nd[rep(seq_len(nrow(nd)), each = length(grid)), ]; long$time <- rep(grid, nrow(nd))
  S <- matrix(predict(fit, newdata = long, type = "surv"), ncol = length(grid), byrow = TRUE)
  10 + rowSums((S[, -1] + S[, -ncol(S)]) / 2) * diff(grid)[1]     # restricted mean age at first birth
}
# mediator: "sec" (binary) or "joint" (sec x union before 18, 4 levels)
mediate <- function(d, med = "sec", zone_mod = FALSE) {
  d$w <- d$wt / mean(d$wt)
  zt <- if (zone_mod) " * zone" else " + zone"
  # Survey round and urban residence are included in both nuisance models.
  adj <- " + religion + urban + factor(year)"
  if (med == "sec") { d$M <- factor(d$sec); mf <- paste0("M ~ A", zt, adj); yf <- paste0("Surv(time, event) ~ A * M", zt, adj) }
  else { d$M <- factor(paste0(d$sec, d$union18)); mf <- paste0("M ~ A", zt, adj)
         yf <- paste0("Surv(time, event) ~ A * M", zt, adj) }
  mm <- multinom(as.formula(mf), d, weights = w, trace = FALSE)
  fy <- stpm2(as.formula(yf), data = d, weights = w, df = 4)
  cdist <- d |> group_by(zone, religion, urban, year) |> summarise(wc = sum(w), .groups = "drop")
  lev <- levels(d$M)
  pm <- function(a) { p <- predict(mm, mutate(cdist, A = a), type = "probs")
                      if (length(lev) == 2) p <- cbind(1 - p, p); colnames(p) <- lev; p }
  Rm <- sapply(lev, function(m) sapply(0:1, function(a)
    list(rmst_profiles(fy, mutate(cdist, A = a, M = factor(m, lev)))))) # [a+1, m] -> vector over cdist
  Q <- function(a, ap, by = NULL) {
    P <- pm(ap); v <- Reduce(`+`, lapply(lev, function(m) P[, m] * Rm[[a + 1, m]]))
    if (is.null(by)) weighted.mean(v, cdist$wc) else c(tapply(v * cdist$wc, cdist[[by]], sum) / tapply(cdist$wc, cdist[[by]], sum))
  }
  top <- tail(lev, 1)  # "1" (secondary) or "10"/"11": the all-secondary scenario uses the secondary level(s)
  sec_lv <- if (med == "sec") "1" else c("10", "11")
  all_sec <- { P <- pm(1); P2 <- P; P2[, ] <- 0
    if (med == "sec") P2[, "1"] <- 1 else { P2[, "10"] <- P[, "00"] + P[, "10"]; P2[, "11"] <- P[, "01"] + P[, "11"] }
    v <- Reduce(`+`, lapply(lev, function(m) P2[, m] * Rm[[2, m]])); weighted.mean(v, cdist$wc) }
  q11 <- Q(1, 1); q10 <- Q(1, 0); q00 <- Q(0, 0)
  r <- tibble(mediator = med, q00 = q00, q10 = q10, q11 = q11, total = q11 - q00, iie = q11 - q10,
              ide = q10 - q00, prop_mediated = (q11 - q10) / (q11 - q00), all_secondary_1990s = all_sec - q11)
  if (zone_mod) {
    z11 <- Q(1, 1, "zone"); z10 <- Q(1, 0, "zone"); z00 <- Q(0, 0, "zone")
    r <- tibble(mediator = med, zone = names(z11), total = z11 - z00, iie = z11 - z10, prop_mediated = (z11 - z10) / (z11 - z00))
  }
  attr(r, "models") <- list(mediator = mm, outcome = fy); r
}
log_step("Point estimates")
pt_sec   <- mediate(d2, "sec")
pt_joint <- mediate(d2, "joint")
pt_zone  <- mediate(d2, "sec", zone_mod = TRUE)
saveRDS(attr(pt_sec, "models"), file.path(o, "sap3_models_primary.rds"))
wcsv(tidy(attr(pt_sec, "models")$mediator), "sap3_mediator_model_coefficients.csv")
wcsv(as_tibble(summary(attr(pt_sec, "models")$outcome)@coef, rownames = "term"), "sap3_outcome_model_coefficients.csv")

# ---- Cluster bootstrap ----
log_step("Bootstrap B = ", B_SAP3)
psus <- unique(d2[, c("psu", "strata")])
bs <- mclapply(seq_len(B_SAP3), function(b) tryCatch({
  mu <- setNames(boot_mult(psus$psu, psus$strata), psus$psu)
  db <- d2[rep(seq_len(nrow(d2)), mu[as.character(d2$psu)]), ]
  st <- function(x) { attr(x, "models") <- NULL; x }
  list(sec = st(mediate(db, "sec")), joint = st(mediate(db, "joint")), zone = st(mediate(db, "sec", TRUE)))
}, error = function(e) structure(conditionMessage(e), class = "try-error")), mc.cores = CORES, mc.set.seed = TRUE)
bs <- bs[!vapply(bs, inherits, logical(1), "try-error")]
if (length(bs) < ceiling(.95 * B_SAP3)) stop("Fewer than 95% of SAP3 bootstrap replicates completed.")
saveRDS(bs, file.path(o, "sap3_bootstrap_replicates.rds"))
ci_tab <- function(pt, reps, by = NULL) {
  r <- bind_rows(reps) |> select(-any_of(c("q00", "q10", "q11"))) |>
    pivot_longer(-any_of(c("mediator", by)), names_to = "estimand") |>
    group_by(across(any_of(c("mediator", by, "estimand")))) |>
    summarise(lo = quantile(value, .025, na.rm = TRUE), hi = quantile(value, .975, na.rm = TRUE), B = n(), .groups = "drop")
  pt |> select(-any_of(c("q00", "q10", "q11"))) |> pivot_longer(-any_of(c("mediator", by)), names_to = "estimand", values_to = "estimate") |>
    left_join(r, by = c("mediator", by, "estimand"))
}
res <- bind_rows(ci_tab(pt_sec, map(bs, "sec")), ci_tab(pt_joint, map(bs, "joint")))
wcsv(res, "sap3_table2_mediation_results.csv")
wcsv(filter(res, mediator == "sec", estimand == "prop_mediated"), "sap3_primary_estimand.csv")
wcsv(ci_tab(pt_zone, map(bs, "zone"), "zone"), "sap3_table3_mediation_by_zone.csv")
wcsv(tibble(cohort = c("1950-1959", "1990-1999"), rmst_age_first_birth = c(pt_sec$q00, pt_sec$q11)), "sap3_rmst_by_cohort.csv")

# ---- Sensitivity: alternative mediator (v106 level only), surveys 2008+ only ----
sens <- bind_rows(
  mediate(mutate(d2, sec = as.integer(medu >= 2)) |> filter(!is.na(sec)), "sec") |> mutate(analysis = "Mediator = v106 secondary+"),
  { x <- filter(d2, year >= 2008 | cohort == "1950-1959"); mediate(x, "sec") |> mutate(analysis = "1990s cohort from 2008+ surveys") },
  mediate(filter(d2, year != 1990), "sec") |> mutate(analysis = "Excluding 1990 survey"))
wcsv(sens, "sap3_table4_sensitivity.csv")

# ---- Figures ----
save_fig(ggplot(filter(cum, !is.na(pct)), aes(cohort, pct, colour = factor(age), group = factor(age))) +
  geom_line(linewidth = 1) + geom_pointrange(aes(ymin = lo, ymax = hi), size = .3) + scale_colour_manual(values = pal) +
  labs(title = "Cumulative probability of first birth by birth cohort", x = "Birth cohort",
       y = "% with a first birth", colour = "By age") + theme_sap, file.path(o, "fig1_cumulative_first_birth.png"), 9, 5)
fd <- filter(res, estimand %in% c("total", "ide", "iie"))
save_fig(ggplot(fd, aes(estimand, estimate, fill = mediator)) + geom_col(position = position_dodge(.8), width = .7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8), width = .25) +
  scale_x_discrete(labels = c(ide = "Direct", iie = "Indirect (via schooling)", total = "Total")) +
  scale_fill_manual(values = pal[c(4, 2)], labels = c(joint = "Schooling + union <18", sec = "Schooling")) +
  labs(title = "Change in restricted mean age at first birth, 1950s vs 1990s cohorts",
       x = NULL, y = "Years (95% CI)", fill = "Mediator") + theme_sap, file.path(o, "fig2_mediation.png"))
log_step("SAP3 done. Proportion mediated = ", round(100 * pt_sec$prop_mediated, 1), "%")