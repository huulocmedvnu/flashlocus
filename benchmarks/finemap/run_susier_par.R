# Throughput: susieR parallel over loci with forked workers (parallel::mclapply), wall time including I/O.
if (nzchar(Sys.getenv("SUSIE_LIB"))) .libPaths(c(Sys.getenv("SUSIE_LIB"), .libPaths()))
suppressMessages({library(susieR); library(parallel)})
a <- commandArgs(TRUE); man <- read.delim(a[1], header = FALSE); cores <- as.integer(a[2])
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(man)), function(i) {
  z <- scan(man[i, 2], quiet = TRUE); p <- length(z)
  R <- matrix(readBin(man[i, 3], "double", p * p, endian = "little"), p, p, byrow = TRUE)
  f <- tryCatch(susie_rss(z = z, R = R, n = man[i, 4], L = 10), error = function(e) NULL)
  if (is.null(f)) NA else length(f$sets$cs)
}, mc.cores = cores, mc.preschedule = FALSE)
cat("susieR", as.character(packageVersion("susieR")), "loci", nrow(man), "cores", cores, "wall_s", as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")
