#!/usr/bin/env Rscript
# dada2_its2.R — DADA2 for Frey ITS2 (PRJNA704912), 234 samples
# ITS-specific: NO truncLen (length-variable amplicon); quality-based filtering only.

suppressMessages(library(dada2))

path <- "trimmed"
out  <- "dada2_out"
dir.create(out, showWarnings = FALSE)

fnFs <- sort(list.files(path, pattern="_1.fastq.gz$", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz$", full.names=TRUE))
sample.names <- sub("_1.fastq.gz$", "", basename(fnFs))
stopifnot(length(fnFs) == 234, length(fnFs) == length(fnRs))

filtFs <- file.path(out, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(out, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- sample.names

ft <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    truncLen=0, maxN=0, maxEE=c(2,2), truncQ=2,
                    minLen=50, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
saveRDS(ft, file.path(out, "filter_stats.rds"))
print(head(ft))

keep   <- file.exists(filtFs)
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF, file.path(out, "errF.rds"))
saveRDS(errR, file.path(out, "errR.rds"))

dadaFs  <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs  <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab    <- makeSequenceTable(mergers)
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus",
                                multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nc, file.path(out, "seqtab_nochim.rds"))

seqs <- getSequences(seqtab.nc)
ids  <- paste0("ASV", seq_along(seqs))
writeLines(paste0(">", ids, "\n", seqs), file.path(out, "asvs.fasta"))

cat(sprintf("\nDONE: %d samples, %d ASVs (post-chimera)\n",
            nrow(seqtab.nc), ncol(seqtab.nc)))
cat(sprintf("ASV length range: %d-%d bp\n",
            min(nchar(seqs)), max(nchar(seqs))))
