#!/usr/bin/env Rscript
# =============================================================================
# run_full_v10.R  (ReactionNormConstraints v10.0.0)
#
# Driver for the full v10 reruns. Each part runs in its own Rscript process with
# a log file, so a failure in one part does not lose the others.
#
#   sim_v9   Experiments 1-25 (the v9.0.0 set), summaries, bootstrap, figures
#            -> ./sim_output_v10/
#   sim_v10  Experiments 26-29 (slope-only QTL class, core simulation,
#            within-environment variance, per-locus variance explained)
#            -> ./sim_output_v10/
#   audit    reaction_norm_audit.R stages 1-4 -> ./audit_output_v10/
#   check    tools/derivation_check.R at full scale
#            -> ./sim_output_v10/derivation_check/
#
# Usage, from the repository root:
#   Rscript run_full_v10.R                         # all parts
#   Rscript run_full_v10.R --parts=sim_v10,audit   # selected parts
#   Rscript run_full_v10.R --smoke                 # SMOKE_TEST mode (~2% replicates)
#   Rscript run_full_v10.R --cores=24              # number of parallel workers
#   Rscript run_full_v10.R --dry-run               # print the plan and estimates only
#
# When sim_v9 and sim_v10 are both requested they run in one process, so
# sim_output_v10/run_metadata_and_session.txt describes the whole simulation run.
#
# Runtime estimates (wall clock, full scale) are in RUNTIME_ESTIMATES below.
# They scale roughly with 1 / (number of workers).
# =============================================================================

# Measured with tools/benchmark_runtime.R (R 4.6.1, reference BLAS, 27 workers on a
# 32-thread Windows workstation): per-replicate cost under full load times the
# replicate count of each experiment, using run_experiment()'s chunking rule.
# The values below are the model's values scaled by 0.7 to match the observed
# Experiment 1 chunk time (12.8 s per 50-task chunk). Treat them as +/- 30%. Re-run
# tools/benchmark_runtime.R for your own machine.
RUNTIME_ESTIMATES <- data.frame(
  part        = c("sim_v9", "sim_v10", "audit", "check"),
  hours_27w   = c(2.1, 0.10, 0.09, 0.05),
  stringsAsFactors = FALSE)
REFERENCE_WORKERS <- 27L

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  if (!grepl("=", hit[1])) return(TRUE)
  sub(paste0("^--", key, "="), "", hit[1])
}
parts   <- strsplit(get_arg("parts", "sim_v9,sim_v10,audit,check"), ",")[[1]]
parts   <- trimws(parts)
smoke   <- isTRUE(get_arg("smoke", FALSE))
dry     <- isTRUE(get_arg("dry-run", FALSE))
cores   <- get_arg("cores", "")
known   <- c("sim_v9", "sim_v10", "audit", "check")
if (any(!parts %in% known)) stop("Unknown part(s): ", paste(setdiff(parts, known), collapse = ", "),
                                 ". Choose from: ", paste(known, collapse = ", "))

# Run from the directory that holds this script.
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!is.na(self) && nzchar(self)) setwd(dirname(normalizePath(self, winslash = "/")))
stopifnot(file.exists("sim_reaction_norm_gwas_v10.R"), file.exists("reaction_norm_audit.R"),
          file.exists(file.path("tools", "derivation_check.R")))

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
workers <- if (nzchar(cores)) as.integer(cores) else max(1L, floor(parallel::detectCores(logical = TRUE) * 0.85))

est <- RUNTIME_ESTIMATES
est$hours_here <- est$hours_27w * REFERENCE_WORKERS / workers
cat("\n  ReactionNormConstraints v10.0.0 full-run driver\n")
cat(sprintf("  Rscript: %s\n  Workers: %d   SMOKE_TEST: %s\n  Parts:   %s\n\n", rscript, workers, smoke,
            paste(parts, collapse = ", ")))
if (!smoke) {
  sel <- est[est$part %in% parts, ]
  if (all(!is.na(sel$hours_here))) {
    cat("  Estimated wall-clock time at full scale:\n")
    for (i in seq_len(nrow(sel))) cat(sprintf("    %-8s %5.1f h\n", sel$part[i], sel$hours_here[i]))
    cat(sprintf("    %-8s %5.1f h\n\n", "total", sum(sel$hours_here)))
  }
}
if (dry) quit(save = "no", status = 0)

dir.create("logs", showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_env <- c(SMOKE_TEST = if (smoke) "TRUE" else "FALSE")
if (nzchar(cores)) base_env <- c(base_env, RN_N_CORES = as.character(workers))

run_part <- function(label, script, env) {
  log <- file.path("logs", sprintf("%s_%s.log", label, stamp))
  cat(sprintf("  [%s] %s started; log: %s\n", format(Sys.time(), "%H:%M:%S"), label, log))
  old <- Sys.getenv(names(env), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    for (nm in names(old)) if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
  }, add = TRUE)
  t0 <- proc.time()[3]
  status <- system2(rscript, script, stdout = log, stderr = log)
  sec <- proc.time()[3] - t0
  cat(sprintf("  [%s] %s finished: exit %d, %.1f min\n", format(Sys.time(), "%H:%M:%S"), label, status, sec / 60))
  data.frame(part = label, exit_status = status, minutes = round(sec / 60, 2), log = log,
             smoke = smoke, workers = workers, stringsAsFactors = FALSE)
}

summary <- list()
sim_parts <- intersect(c("sim_v9", "sim_v10"), parts)
if (length(sim_parts)) {
  env <- c(base_env, RN_RUN_V9_SET = if ("sim_v9" %in% parts) "TRUE" else "FALSE",
           RN_RUN_V10_SET = if ("sim_v10" %in% parts) "TRUE" else "FALSE")
  summary[[length(summary) + 1]] <- run_part(paste(sim_parts, collapse = "+"), "sim_reaction_norm_gwas_v10.R", env)
}
if ("audit" %in% parts)
  summary[[length(summary) + 1]] <- run_part("audit", "reaction_norm_audit.R", base_env)
if ("check" %in% parts) {
  env <- c(base_env, RN_CHECK_OUT = file.path("sim_output_v10", "derivation_check"),
           RN_CHECK_REPS_A = if (smoke) "4" else "200",
           RN_CHECK_REPS_B = if (smoke) "6" else "300",
           RN_CHECK_REPS_C = if (smoke) "4" else "200")
  dir.create("sim_output_v10", showWarnings = FALSE)
  summary[[length(summary) + 1]] <- run_part("check", file.path("tools", "derivation_check.R"), env)
}
summary <- do.call(rbind, summary)
sf <- file.path("logs", sprintf("run_summary_%s.csv", stamp))
write.csv(summary, sf, row.names = FALSE)
cat("\n"); print(summary[, c("part", "exit_status", "minutes")], row.names = FALSE)
cat(sprintf("\n  Summary written to %s\n", sf))
if (any(summary$exit_status != 0)) quit(save = "no", status = 1)
