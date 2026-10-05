# Shared helpers for the SoftwareX validation and benchmark scripts.
# Every script writes a CSV to results/ and a sidecar .meta.txt that records
# the package versions, the hexify commit and the run time.

suppressPackageStartupMessages({
  library(hexify)
  library(sf)
})

results_dir <- function() {
  d <- file.path(dirname(sys_script_path()), "results")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

sys_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) normalizePath(f) else normalizePath("paper/bench/bench_common.R")
}

# Points distributed uniformly over the sphere.
sphere_points <- function(n, seed) {
  set.seed(seed)
  data.frame(lon = runif(n, -180, 180),
             lat = asin(runif(n, -1, 1)) * 180 / pi)
}

# Great-circle distance in km on the sphere of the authalic Earth radius.
gc_km <- function(lon1, lat1, lon2, lat2, radius = 6371.007180918475) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * radius * asin(pmin(1, sqrt(a)))
}

cpu_name <- function() {
  out <- tryCatch(switch(Sys.info()[["sysname"]],
    Windows = system2("powershell", c("-NoProfile", "-Command",
                                      shQuote("(Get-CimInstance Win32_Processor).Name")),
                      stdout = TRUE),
    Darwin = system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE),
    sub("^.*: ", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1])),
    error = function(e) NA_character_)
  trimws(out[1])
}

write_result <- function(df, name, extra = character()) {
  out <- file.path(results_dir(), paste0(name, ".csv"))
  write.csv(df, out, row.names = FALSE)
  pkgs <- c("hexify", "dggridR", "h3o", "h3r", "sf", "s2")
  have <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  commit <- tryCatch(
    system2("git", c("-C", shQuote(dirname(sys_script_path())), "rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_)
  dirty <- tryCatch(
    length(system2("git", c("-C", shQuote(dirname(sys_script_path())), "status", "--porcelain",
                            "--", "../../R", "../../src", "."), stdout = TRUE, stderr = FALSE)) > 0,
    error = function(e) NA)
  meta <- c(
    paste("script:", basename(sys_script_path())),
    paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("hexify commit:", commit),
    paste("uncommitted changes in R/, src/ or paper/bench/:", dirty),
    paste("hexify installed version:", as.character(packageVersion("hexify"))),
    paste("R:", R.version.string),
    paste("platform:", R.version$platform),
    paste("cpu:", cpu_name()),
    paste0(have, ": ", vapply(have, function(p) as.character(packageVersion(p)), "")),
    extra)
  writeLines(meta, sub("\\.csv$", ".meta.txt", out))
  message("wrote ", out)
  invisible(out)
}
