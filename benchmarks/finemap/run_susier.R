# Reference: susieR::susie_rss(z, R, n, L = 10) with defaults; time excludes file I/O.
if (!grepl("renv", R.home())) .libPaths(c("${BENCH_ROOT:-.}/Rlib", .libPaths()))
if (nzchar(Sys.getenv("SUSIE_LIB"))) .libPaths(c(Sys.getenv("SUSIE_LIB"), .libPaths()))
cat("BLAS:", sessionInfo()$BLAS, "\n")
suppressMessages(library(susieR))
cat("susieR", as.character(packageVersion("susieR")), "Rfast", requireNamespace("Rfast", quietly = TRUE), "\n")
a <- commandArgs(TRUE); man <- read.delim(a[1], header = FALSE); out <- a[2]; dir.create(out, showWarnings = FALSE)
tim <- c()
for (i in seq_len(nrow(man))) {
  id <- man[i, 1]; z <- scan(man[i, 2], quiet = TRUE); p <- length(z)
  R <- matrix(readBin(man[i, 3], "double", p * p, endian = "little"), p, p, byrow = TRUE); n <- man[i, 4]
  t <- system.time(f <- tryCatch(susie_rss(z = z, R = R, n = n, L = 10), error = function(e) e))[["elapsed"]]
  if (inherits(f, "error")) { tim <- rbind(tim, data.frame(id = id, p = p, seconds = NA, niter = NA, n_cs = NA, err = conditionMessage(f))); next }
  writeLines(sprintf("%.10e", f$pip), file.path(out, paste0(id, ".pip")))
  cs <- f$sets$cs
  if (length(cs)) writeLines(paste(sub("L", "", names(cs)), sapply(cs, function(x) paste(x, collapse = ",")), sep = "\t"), file.path(out, paste0(id, ".cs")))
  else file.create(file.path(out, paste0(id, ".cs")))
  write.table(f$lbf_variable, file.path(out, paste0(id, ".lbf_variable.tsv")), sep = "\t", row.names = FALSE, col.names = FALSE)
  tim <- rbind(tim, data.frame(id = id, p = p, seconds = t, niter = f$niter, n_cs = length(cs), err = ""))
}
write.table(tim, file.path(out, "timing.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
