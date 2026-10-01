#!/usr/bin/env Rscript
# Exact comparison of a nanorust alldata.rds with an R parser output (any layout:
# v18_Br/E, v18_BrE, DS421_v2). Columns only present in the R file at step 03
# (signal, mapping_table) are dropped. Rows are ordered by (chrom, read_id, flag, start).
# Usage: compare_any.R new.rds ref.rds
suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
new <- readRDS(args[1])
ref <- readRDS(args[2]) %>% select(-any_of(c("signal", "mapping_table")))
ok <- TRUE
report <- function(what, pass, detail = "") {
  cat(sprintf("%-28s %s %s\n", what, if (pass) "OK  " else "FAIL", detail))
  if (!pass) ok <<- FALSE
}
report("columns", identical(names(new), names(ref)),
       paste0("new: ", paste(names(new), collapse = ","), if (!identical(names(new), names(ref))) paste0(" | ref: ", paste(names(ref), collapse = ",")) else ""))
report("nrow", nrow(new) == nrow(ref), sprintf("%d vs %d", nrow(new), nrow(ref)))
o <- function(d) d %>% arrange(chrom, read_id, flag, start)
new <- o(new); ref <- o(ref)
key <- function(d) paste(d$read_id, d$flag, d$chrom, d$start)
only_new <- setdiff(key(new), key(ref)); only_ref <- setdiff(key(ref), key(new))
report("same mappings", !length(only_new) && !length(only_ref), sprintf("only_new=%d only_ref=%d", length(only_new), length(only_ref)))
if (length(only_new)) print(head(new[key(new) %in% only_new, 1:6]))
if (length(only_ref)) print(head(ref[key(ref) %in% only_ref, 1:6]))
if (identical(key(new), key(ref))) {
  n <- new; r <- ref  # same rows in the same (stable) order: compare row by row
} else {
  common <- intersect(key(new), key(ref))
  n <- new[match(common, key(new)), ]; r <- ref[match(common, key(ref)), ]
}
for (col in intersect(names(n), names(r))) {
  if (col == "signalbin") next
  same <- identical(n[[col]], r[[col]])
  detail <- ""
  if (!same && is.numeric(n[[col]])) {
    bad <- which(!(n[[col]] == r[[col]] | (is.na(n[[col]]) & is.na(r[[col]]))) | is.na(n[[col]]) != is.na(r[[col]]))
    detail <- sprintf("%d differ, max |diff| %.3g", length(bad), max(c(0, abs(n[[col]] - r[[col]])), na.rm = TRUE))
  }
  report(col, same, detail)
}
sb_same <- mapply(identical, n$signalbin, r$signalbin)
report("signalbin (nested, exact)", all(sb_same), sprintf("%d of %d differ", sum(!sb_same), length(sb_same)))
if (!all(sb_same)) {
  i <- which(!sb_same)[1]
  cat("first differing signalbin, row", i, n$read_id[i], "\n")
  print(all.equal(n$signalbin[[i]], r$signalbin[[i]]))
}
cat(if (ok) "ALL IDENTICAL\n" else "DIFFERENCES FOUND\n")
quit(status = if (ok) 0 else 1)
