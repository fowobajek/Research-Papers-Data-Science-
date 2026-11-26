# Shared configuration for the SAP analysis pipeline (SAP1-SAP3)
# Reuses the EDA helpers (survey design, file finder) and the harmonisation rules
# (HW anthropometry links, 1990 WI wealth link, 1990 zones from SSTATE).
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(survey); library(data.table); library(broom)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, warn = 1)

SEED   <- 20261020L                                         # pre-registered seed (all SAPs)
B_SAP2 <- as.integer(Sys.getenv("BOOT_B_SAP2", "1000"))     # cluster bootstrap replicates, SAP2
B_SAP3 <- as.integer(Sys.getenv("BOOT_B_SAP3", "1000"))     # cluster bootstrap replicates, SAP3
MI_M   <- as.integer(Sys.getenv("MI_M", "20"))              # multiple imputations
CORES  <- as.integer(Sys.getenv("CORES", "3"))

ROOT <- normalizePath(Sys.getenv("SAP_ROOT", "."), mustWork = FALSE)
PROJECT_ROOT <- normalizePath(file.path(ROOT, ".."), mustWork = FALSE)
DATA_DIR <- Sys.getenv("DHS_DIR", file.path(PROJECT_ROOT, "data", "unpacked"))
# Allow a complete, immutable release to be written alongside the working
# outputs without overwriting them.  The default preserves the established
# interactive workflow; release builds set SAP_OUT explicitly.
OUT <- Sys.getenv("SAP_OUT", file.path(ROOT, "output"))
for (s in c("data", "sap1", "sap2", "sap3", "reproducibility"))
  dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)

num <- function(x) as.numeric(zap_labels(x))
lab <- function(x) as.character(as_factor(x))
list_dhs <- function() {
  f <- list.files(DATA_DIR, pattern = "\\.sav$", recursive = TRUE, ignore.case = TRUE, full.names = TRUE)
  if (!length(f)) stop("No extracted .sav files found in ", normalizePath(DATA_DIR, mustWork = FALSE),
                       ". Run Rscript scripts/00_materialize_data.R from the project root, or set DHS_DIR.")
  tibble(path = f, year = as.integer(basename(dirname(dirname(f)))),
         recode = substr(toupper(basename(f)), 3, 4))
}
find_file <- function(yr, rc) {
  x <- filter(list_dhs(), year == yr, recode == rc)
  if (nrow(x) > 1) stop("Expected one ", rc, " file for ", yr, "; found ", nrow(x))
  if (nrow(x)) x$path[1] else NA_character_
}
read_vars <- function(path, vars) {
  d <- read_sav(path, col_select = any_of(c(vars, toupper(vars)))); names(d) <- tolower(names(d))
  for (v in setdiff(vars, names(d))) d[[v]] <- NA_real_
  d
}
YEARS  <- c(1990L, 2003L, 2008L, 2013L, 2018L, 2024L)
WEALTH <- c("Poorest", "Poorer", "Middle", "Richer", "Richest")
ZONES  <- c("North Central", "North East", "North West", "South East", "South South", "South West")
ZONE_1990 <- c(Abuja = "North Central", Benue = "North Central", Kwara = "North Central",
  Niger = "North Central", Plateau = "North Central", Bauchi = "North East", Borno = "North East",
  Gongola = "North East", Kaduna = "North West", Kano = "North West", Katsina = "North West",
  Sokoto = "North West", Anambra = "South East", Imo = "South East", `Akwa Ibom` = "South South",
  Bendel = "South South", `Cross River` = "South South", Rivers = "South South",
  Lagos = "South West", Ogun = "South West", Ondo = "South West", Oyo = "South West")

to_zone <- function(d, yr) {
  if (yr == 1990) {
    z <- unname(ZONE_1990[trimws(lab(d$sstate))])
    if (anyNA(z)) stop("Unmapped 1990 state")
    return(z)
  }
  z <- tools::toTitleCase(tolower(lab(d$v024))); z[!z %in% ZONES] <- NA; z
}
# 1990 wealth quintile from WI file: WHHID = CASEID[1:12]
merge_wi <- function(d, yr) {
  f <- find_file(yr, "WI")
  if (is.na(f) || !all(is.na(d$v190))) return(d)
  wi <- read_sav(f); names(wi) <- tolower(names(wi))
  wi <- transmute(wi, .hh = whhid, .q = num(wlthind5), .s = num(wlthindf))
  out <- left_join(mutate(d, .hh = substr(caseid, 1, 12)), wi, by = ".hh")
  stopifnot(nrow(out) == nrow(d))
  out$v190 <- out$.q; out$wealth_score <- out$.s; select(out, -.hh, -.q, -.s)
}
# 1990/2003 anthropometry from HW files (DHS MERGE.DOC keys)
merge_hw <- function(d, yr) {
  f <- find_file(yr, "HW")
  # Retain whether an early-round KR record found a linked HW record.  This
  # distinguishes a failed link from a linked child with no usable HAZ.
  if (is.na(f) || !all(is.na(num(d$hw70)))) {
    d$hw_file_match <- NA
    return(d)
  }
  hw <- read_sav(f); names(hw) <- tolower(names(hw))
  hw <- transmute(hw, across(any_of(c("hwcaseid", "hwhhid"))), .k2 = num(hwline),
                  .hc = num(hc70), hw_file_match = TRUE)
  d <- select(d, -hw70)
  out <- if ("hwcaseid" %in% names(hw))
    left_join(mutate(d, .k1 = caseid, .k2 = num(midx)), rename(hw, .k1 = hwcaseid), by = c(".k1", ".k2")) else
    left_join(mutate(d, .k1 = substr(caseid, 1, 12), .k2 = num(b16)), rename(hw, .k1 = hwhhid), by = c(".k1", ".k2"))
  stopifnot(nrow(out) == nrow(d))
  rename(out, hw70 = .hc) |> select(-.k1, -.k2)
}
# Common design variables; cluster id unique across surveys
add_design <- function(d, yr) {
  d$wt <- num(d$v005) / 1e6
  d$psu <- yr * 1e5 + num(d$v021)                       # numeric ids, unique across surveys
  d$strata <- if (length(unique(na.omit(num(d$v022)))) > 1) yr * 1e3 + num(d$v022) else
    yr * 1e3 + num(d$v024) * 10 + num(d$v025)
  d
}
# Weighted fractional rank (ridit) of wealth within survey, ties at midpoint
ridit <- function(q, w) {
  ok <- !is.na(q); r <- rep(NA_real_, length(q))
  tw <- tapply(w[ok], q[ok], sum); cw <- cumsum(tw) / sum(tw); mid <- cw - (tw / sum(tw)) / 2
  r[ok] <- mid[as.character(q[ok])]; r
}
# Submission figures need print-quality source files.  PNGs remain convenient
# inputs for DOCX while the release packager converts the final main figures to
# TIFF/PDF as required by the journal.
save_fig <- function(p, path, w = 8, h = 5) ggsave(path, p, width = w, height = h, dpi = 300, bg = "white")
theme_sap <- theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"),
                                                   panel.grid.minor = element_blank())
pal <- c("#1b5e20", "#c49a00", "#8d2b0b", "#2e6f95", "#6b6b6b", "#5a3e85", "#b5651d")
# Cluster bootstrap: resample PSUs within strata; returns per-row multiplicity
boot_mult <- function(psu, strata) {
  u <- unique(data.frame(psu, strata))
  m <- unlist(lapply(split(u$psu, u$strata), function(p) table(factor(sample(p, length(p), TRUE), levels = p))))
  names(m) <- sub("^[^.]*\\.", "", names(m))
  as.numeric(m[as.character(psu)])
}
log_step <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
cf <- function(x) as.numeric(coef(x)); se_ <- function(x) as.numeric(SE(x))

# Survey-matched boundaries are used only for post-estimation mapping.  The 1990
# archive has four non-comparable regions, so it is deliberately excluded from
# six-zone maps.  Neither the 1999 nor the 2025 archive belongs to this project.
load_zone_boundaries <- function(yr) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for spatial outputs.")
  if (!yr %in% c(2003L, 2008L, 2013L, 2018L, 2024L))
    stop("Six-zone boundaries are only comparable for 2003-2024; do not map 1990 with later zones.")
  z <- file.path(PROJECT_ROOT, "shape files",
                 sprintf("Nigeria %d DHS - sdr_subnational_boundaries_2026-10-22.zip", yr))
  if (!file.exists(z)) stop("Missing survey-matched boundary archive: ", z)
  d <- file.path(OUT, "cache", "boundaries", as.character(yr)); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  shp <- file.path(d, "shps", "sdr_subnational_boundaries.shp")
  if (!file.exists(shp)) utils::unzip(z, exdir = d)
  x <- sf::st_read(shp, quiet = TRUE)
  if (!all(x$SVYYEAR == yr) || !all(x$REGVAR == "hv024")) stop("Boundary metadata does not match DHS ", yr)
  if (any(!sf::st_is_valid(x, NA_on_exception = TRUE))) x <- sf::st_make_valid(x)
  x |> dplyr::transmute(zone = tools::toTitleCase(tolower(REGNAME)), geometry = geometry)
}

write_zone_map <- function(results, yr, value, path, title) {
  b <- load_zone_boundaries(yr)
  d <- dplyr::left_join(b, results |> dplyr::mutate(zone = tools::toTitleCase(tolower(zone))), by = "zone")
  if (all(is.na(d[[value]]))) stop("No estimates joined to boundaries for ", yr)
  p <- ggplot(d) + geom_sf(aes(fill = .data[[value]]), colour = "white", linewidth = .2) +
    scale_fill_viridis_c(na.value = "grey90") + labs(title = title, fill = value) + theme_void()
  save_fig(p, path, 7, 6)
}
if (Sys.getenv("OPENBLAS_NUM_THREADS") == "") warning("Start R via run.sh (OPENBLAS_NUM_THREADS=1) to avoid CPU oversubscription in bootstrap workers")