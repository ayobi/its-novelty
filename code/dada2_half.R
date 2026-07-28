#!/usr/bin/env Rscript
suppressMessages(library(dada2))
args <- commandArgs(trailingOnly=TRUE)
samplelist <- args[1]; out <- args[2]
dir.create(out, showWarnings=FALSE, recursive=TRUE)

sn <- readLines(samplelist)
fnFs <- file.path("trimmed", paste0(sn, "_1.fastq.gz"))
fnRs <- file.path("trimmed", paste0(sn, "_2.fastq.gz"))
stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))

filtFs <- file.path(out,"filtered",paste0(sn,"_F.fastq.gz"))
filtRs <- file.path(out,"filtered",paste0(sn,"_R.fastq.gz"))
names(filtFs) <- names(filtRs) <- sn

filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=0, maxN=0,
              maxEE=c(2,2), truncQ=2, minLen=50, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
keep <- file.exists(filtFs); filtFs<-filtFs[keep]; filtRs<-filtRs[keep]

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
seqtab <- makeSequenceTable(mergers)
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab.nc, file.path(out,"seqtab_nochim.rds"))

seqs <- getSequences(seqtab.nc)
writeLines(paste0(">ASV", seq_along(seqs), "\n", seqs), file.path(out,"asvs.fasta"))
cat(sprintf("DONE %s: %d samples, %d ASVs, len %d-%d\n",
    out, nrow(seqtab.nc), ncol(seqtab.nc), min(nchar(seqs)), max(nchar(seqs))))
