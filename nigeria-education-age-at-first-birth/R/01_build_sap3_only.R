# SAP3-only harmonisation of DHS Individual Recode data for age at first birth
source(file.path(Sys.getenv("SAP_ROOT", "."), "R", "00_config.R"))
fix_year <- function(y) ifelse(y < 100, 1900 + y, y)
religion <- function(x) {
  l <- tolower(lab(x))
  case_when(grepl("islam|muslim|moslem", l) ~ "Muslim",
            grepl("catholic|protestant|christian|pentecost|evangel|anglican|methodist", l) ~ "Christian",
            is.na(x) ~ NA_character_, TRUE ~ "Other/none")
}
counts <- list()

# ---------------- SAP3: women, age at first birth ----------------
ir <- map_dfr(YEARS, function(yr) {
  d <- read_vars(find_file(yr, "IR"), c("caseid","v001","v005","v008","v010","v011","v012","v021","v022",
               "v024","v025","v106","v133","v130","v190","v201","v211","v212","v509","v511","sstate"))
  d <- merge_wi(d, yr) |> add_design(yr)
  agefb <- ifelse(num(d$v201) > 0, (num(d$v211) - num(d$v011)) / 12, NA)
  tibble(year = yr, psu = d$psu, strata = d$strata, wt = d$wt, byear = fix_year(num(d$v010)),
         age = (num(d$v008) - num(d$v011)) / 12, event = as.integer(num(d$v201) > 0),
         t_fb = ifelse(num(d$v201) > 0, agefb, (num(d$v008) - num(d$v011)) / 12),
         edyrs = ifelse(num(d$v133) > 30, NA, num(d$v133)), medu = num(d$v106),
         afu = ifelse(num(d$v511) > 0 & num(d$v511) < 98, num(d$v511), NA),
         ever_union = !is.na(num(d$v511)) & num(d$v511) < 98,
         urban = as.integer(num(d$v025) == 1), religion = religion(d$v130), zone = to_zone(d, yr))
}) |> mutate(
  cohort = cut(byear, seq(1940, 2010, 10), paste0(seq(1940, 2000, 10), "-", seq(1949, 2009, 10)), right = FALSE),
  # Observed educational-attainment proxy, not a retrospective measure at age 15.
  sec = as.integer(edyrs >= 7 | medu >= 2),
  # Second mediator: first union before 18 (never-in-union = no)
  union18 = as.integer(!is.na(afu) & afu < 18),
  t_fb = pmax(t_fb, 10.01))
ir$sec[is.na(ir$edyrs) & is.na(ir$medu)] <- NA
saveRDS(ir, file.path(OUT, "data", "sap3_women.rds"))
counts$sap3 <- ir |> group_by(year) |> summarise(women = n(), first_births = sum(event),
  missing_school = sum(is.na(sec)), missing_cohort = sum(is.na(cohort)))
log_step("SAP3 data: ", nrow(ir), " women")
iwalk(counts, ~ write_csv(.x, file.path(OUT, "data", paste0(.y, "_sample_counts.csv"))))