#!/usr/bin/env Rscript
# asv_abundance.R -- per-ASV prevalence and read totals from the DADA2 table.
#
#   Rscript asv_abundance.R [seqtab.rds] [asvs.fasta] [out.tsv]
#
# defaults: dada2_out/seqtab_nochim.rds, dada2_out/asvs.fasta, asv_abundance.tsv
#
# Writes asv_id, n_samples, total_reads. Maps FASTA records to seqtab columns by
# sequence identity first, since that is verifiable, and falls back to column order
# only if the counts match exactly. No re-denoising.

args <- commandArgs(TRUE)
st_path  <- if (length(args) >= 1) args[1] else "dada2_out/seqtab_nochim.rds"
fa_path  <- if (length(args) >= 2) args[2] else "dada2_out/asvs.fasta"
out_path <- if (length(args) >= 3) args[3] else "asv_abundance.tsv"

st <- readRDS(st_path)
cat("seqtab:", nrow(st), "samples x", ncol(st), "sequences\n")

lines <- readLines(fa_path)
h <- grep("^>", lines)
if (!length(h)) stop("no FASTA headers in ", fa_path)
ids <- sub("^>", "", sub(" .*", "", lines[h]))
starts <- h + 1
ends <- c(h[-1] - 1, length(lines))
seqs <- vapply(seq_along(ids),
               function(i) if (starts[i] > ends[i]) ""
                           else paste0(lines[starts[i]:ends[i]], collapse = ""),
               character(1))
cat("fasta: ", length(ids), " records, first id ", ids[1], "\n", sep = "")

idx <- match(toupper(seqs), toupper(colnames(st)))
matched <- sum(!is.na(idx))
cat("matched by sequence:", matched, "of", length(ids), "\n")
if (matched < 0.99 * length(ids)) {
  if (length(ids) == ncol(st)) {
    cat("sequences differ (trimmed?), falling back to column order\n")
    idx <- seq_along(ids)
  } else {
    stop("cannot map: ", length(ids), " records against ", ncol(st), " columns. ",
         "Try the other FASTA in dada2_out/.")
  }
}

keep <- !is.na(idx)
sub_st <- st[, idx[keep], drop = FALSE]
d <- data.frame(asv_id      = ids[keep],
                n_samples   = as.integer(colSums(sub_st > 0)),
                total_reads = as.integer(colSums(sub_st)))
write.table(d, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("wrote", out_path, "\n")
cat("  in 1 sample only: ", sum(d$n_samples == 1), "\n", sep = "")
cat("  in >=2 samples:   ", sum(d$n_samples >= 2), "\n", sep = "")
cat("  with >=10 reads:  ", sum(d$total_reads >= 10), "\n", sep = "")
cat("  both filters:     ", sum(d$n_samples >= 2 & d$total_reads >= 10), "\n", sep = "")
