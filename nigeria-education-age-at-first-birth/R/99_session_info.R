# Save exact R environment after running the SAP3 pipeline.
source(file.path(Sys.getenv("SAP_ROOT", "."), "R", "00_config.R"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "reproducibility", "R_session_info.txt"))
