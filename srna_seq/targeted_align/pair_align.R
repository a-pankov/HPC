args <- commandArgs(TRUE)
samp <- args[1]

library(ShortRead)
library(Biostrings)

ref <- readRNAStringSet('/songlab/shared/hairpin.fa.gz')
#ref <- readRNAStringSet('/songlab/shared/mature.fa.gz')
mumu_inds <- grep('Mus musculus', names(ref))

ref_mumu <- ref[mumu_inds]

ref_mumu_dna <- DNAStringSet(ref_mumu)

uniq_refs <- aggregate( names(ref_mumu_dna) ~ as.character(ref_mumu_dna), FUN = paste0, collapse = ';')

uniq_refs_dna <- DNAStringSet(uniq_refs[,1])
names(uniq_refs_dna) <- uniq_refs[,2]

quers <- readFastq(paste0(samp, '_noAdap_noArchAdap3_5_nn_nn.fastq.gz'))

tbl_quers <- table(sread(quers))
quers_uniq <- DNAStringSet(names(tbl_quers))
count_uniq <- as.vector(tbl_quers)

rm(quers, tbl_quers)

library(parallel)
## lscores_align <- mclapply(quers_uniq, function(x) pairwiseAlignment(uniq_refs_dna, x, type = 'global-local', scoreOnly = T), mc.cores = 64)

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -Inf, baseOnly = TRUE)
lscores_align <- mclapply(quers_uniq, function(x) pairwiseAlignment( uniq_refs_dna, x, type = "overlap", scoreOnly= T, substitutionMatrix = mat, gapOpening = -Inf), mc.cores = 64)


save(lscores_align, quers_uniq, count_uniq, file = paste0(samp, "_srna_align.RData"))



