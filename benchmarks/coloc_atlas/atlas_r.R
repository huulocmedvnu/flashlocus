# Reference pipeline: same atlas with coloc 5.2.3 coloc.bf_bf (what coloc.susie calls), data.table fread with all
# threads, mclapply over loci. Usage: Rscript atlas_r.R loci.tsv name lbf.gz cs.gz out.tsv cores
suppressPackageStartupMessages({ library(data.table); library(coloc); library(parallel) })
a <- commandArgs(TRUE); cores <- as.integer(a[6]); setDTthreads(cores)
t0 <- Sys.time()
L <- fread(a[1], header = FALSE, col.names = c("id", "fit", "keys"))
loci <- lapply(seq_len(nrow(L)), function(r) {
  k <- readLines(L$keys[r]); k[k == "NA"] <- paste0("NA#", which(k == "NA") - 1)
  cs <- fread(paste0(L$fit[r], ".cs.tsv")); if (!nrow(cs)) return(NULL)
  lbf <- as.matrix(fread(paste0(L$fit[r], ".lbf_variable.tsv"), header = FALSE))
  bf <- lbf[cs$effect, , drop = FALSE]; colnames(bf) <- k
  ok <- !startsWith(k, "NA#"); if (!any(ok)) return(NULL); s <- tstrsplit(k[ok], "_", fixed = TRUE)
  chr <- s[[1]][1]; pos <- as.numeric(s[[2]][s[[1]] == chr])
  list(id = L$id[r], chr = chr, lo = min(pos), hi = max(pos), bf = bf, eff = cs$effect)
})
loci <- Filter(Negate(is.null), loci)
cs <- fread(a[4], select = c("molecular_trait_id", "cs_id", "region"))
cs[, L := as.integer(sub(".*_L", "", cs_id))]
g <- unique(cs[, .(molecular_trait_id, region)])
g[, `:=`(chr = sub("^chr", "", sub(":.*", "", region)), lo = pmax(0, as.numeric(sub("^.*:(-?[0-9]+)-[0-9]+$", "\\1", region))),
         hi = as.numeric(sub(".*-", "", region)))]
lt <- rbindlist(lapply(loci, function(x) data.table(chr = x$chr, llo = x$lo, lhi = x$hi)))
g <- g[sapply(seq_len(nrow(g)), function(i) any(lt$chr == g$chr[i] & g$lo[i] <= lt$lhi & g$hi[i] >= lt$llo))]
gl <- cs[molecular_trait_id %in% g$molecular_trait_id, .(L = sort(unique(L))), by = molecular_trait_id]
lbf <- fread(a[3], nThread = cores)
lbf <- lbf[molecular_trait_id %in% g$molecular_trait_id]
v <- tstrsplit(lbf$variant, "_", fixed = TRUE)
lbf[, key := paste(sub("^chr", "", v[[1]]), v[[2]], pmin(v[[3]], v[[4]]), pmax(v[[3]], v[[4]]), sep = "_")]
genes <- lapply(split(lbf, by = "molecular_trait_id", keep.by = TRUE), function(d) {
  id <- d$molecular_trait_id[1]; ls <- gl[molecular_trait_id == id, L]
  k <- d$key; dup <- duplicated(k); k[dup] <- paste0("DUP#", which(dup) - 1)
  m <- t(as.matrix(d[, paste0("lbf_variable", ls), with = FALSE])); colnames(m) <- k
  gi <- g[molecular_trait_id == id]
  list(id = id, chr = gi$chr, lo = gi$lo, hi = gi$hi, bf = m, ls = ls)
})
t1 <- Sys.time()
res <- mclapply(loci, function(x) {
  out <- list()
  for (gg in genes) {
    if (gg$chr != x$chr || gg$lo > x$hi || gg$hi < x$lo) next
    r <- suppressWarnings(coloc.bf_bf(x$bf, gg$bf))
    if (is.data.table(r) && "nsnps" %in% names(r) && is.na(r$nsnps[1])) next
    s <- r$summary
    out[[length(out) + 1]] <- data.table(locus = x$id, tissue = a[2], gene = gg$id, gwas_L = x$eff[s$idx1], qtl_L = gg$ls[s$idx2],
      nsnps = s$nsnps, hit1 = s$hit1, hit2 = s$hit2, PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
      PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf)
  }
  rbindlist(out)
}, mc.cores = cores)
res <- rbindlist(res)
fwrite(res, a[5], sep = "\t", na = "NaN")
t2 <- Sys.time()
cat(sprintf("R atlas %s: %d genes, %d pairs, load %.1f s, coloc %.1f s, total %.1f s\n", a[2], length(genes), nrow(res),
            as.numeric(difftime(t1, t0, units = "secs")), as.numeric(difftime(t2, t1, units = "secs")), as.numeric(difftime(t2, t0, units = "secs"))))
