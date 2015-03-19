args <- commandArgs(TRUE)
samp <- args[1]

#setwd("/home/apankov/ilf2_julia/rib_pro/tophat")

#library(BSgenome)
library(GenomicFeatures)
load("gtf_parsed.RData")
#load("trans_seq.RData")
library(parallel)

#dat <- read.table(file = 'JF001_index12_CTTGTA_L005_R1_001_cov.txt', F, sep = '\t', as.is = T)

dat <- read.table(file = paste0('../',samp,"_cov.txt") , F, sep = '\t', as.is = T)
dat_GR <- with(dat, GRanges(seqnames = Rle(V1), ranges = IRanges(V2+1, V3), cov=V4) )
dat_cov <- coverage(dat_GR, weight = dat$V4)

#subset_GR <- subset(exons_GR, seqnames(exons_GR) == "chr13")

exons_GR_split <- split(exons_GR, seqnames(exons_GR))

trans_cov_split <- mclapply(exons_GR_split, function(x) dat_cov[x], mc.cores = 20)
trans_cov <- Reduce(c, trans_cov_split)

#trans_cov <- dat_cov[exons_GR]
save(trans_cov, file = paste0('./',samp,"_exon_cov.RData"))

trans_cov_ordr <- trans_cov[with(gtf_subset, order(transcript_id, as.numeric(exon_number)))]
gtf_subset_ordr <- gtf_subset[with(gtf_subset, order(transcript_id, as.numeric(exon_number))),]

trans_cov_ordr_split <- split(trans_cov_ordr, gtf_subset_ordr$transcript_id )
rev_inds <- which(unique(gtf_subset_ordr[, c("transcript_id", "strand")])$strand == '-')

trans_cov_ordr_split[rev_inds] <- mclapply(trans_cov_ordr_split[rev_inds], function(x) lapply(x, function(y) rev(y) ), mc.cores = 20)

trans_cov_comb <- mclapply(trans_cov_ordr_split, function(x) Reduce(c, x), mc.cores = 20)
save(trans_cov_comb, file = paste0('./',samp,"_trans_cov.RData"))


