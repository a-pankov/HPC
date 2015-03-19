setwd("/home/apankov/ilf2_julia/rib_pro/tophat")
#anno_full <- read.delim("methyl/HumanMethylation450_15017482_v.1.2.csv", header=TRUE, sep=",", as.is=TRUE)

library(GenomicFeatures)
#anno_good <- anno_full[!is.na(anno_full$MAPINFO),]
#anno_good_GR <- with(anno_good, GRanges(seqnames = Rle(paste0("chr", CHR)), ranges = IRanges(MAPINFO, MAPINFO)))


myfile <- "gencode.vM1.annotation_cuffclean_verified.gtf"
library(stringr)
gtf <- read.delim(myfile, header=FALSE, as.is = T)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
chronly <- c(paste0("chr",1:22), "X", "Y")
gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows

#gtf_exon <- subset(gtf, feature == "exon")
gtf_cds <- subset(gtf, feature == "CDS")
#  gtf[, length(gtf[1,])] <- gsub(" {1,}|_id|;", "", gtf$attributes)
#  gtf[, length(gtf[1,])] <- gsub(" {1,}|gene_id|;.*", "", gtf$attributes)

gtf_subset <- gtf_cds

atts <- strsplit(gtf_subset$attributes, ";")
atts_info <- sapply(atts, function(x) gsub("gene_id|transcript_id|exon_id|exon_number", "", grep("gene_id|transcript_id|exon_id|exon_number", x, value = T)))
atts_info_clean <- t(gsub("^\\s+|\\s+$", "", atts_info))
colnames(atts_info_clean) <- c('gene_id', 'transcript_id', 'exon_number', 'exon_id')
gtf_subset <- cbind.data.frame(gtf_subset, atts_info_clean, stringsAsFactors = F)

library(BSgenome.Mmusculus.UCSC.mm9)
exons_GR <- with(gtf_subset, GRanges(seqnames = Rle(seqname), ranges = IRanges(start, end), strand = Rle(strand)) )
save(exons_GR, gtf_subset, file = 'byCDS/gtf_parsed.RData')

exon_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm9, exons_GR)

exon_seqs_ordr <- exon_seqs[with(gtf_subset, order(transcript_id, as.numeric(exon_number)))]
gtf_subset_ordr <- gtf_subset[with(gtf_subset, order(transcript_id, as.numeric(exon_number))),]

exon_seqs_ordr_split <- split(exon_seqs_ordr, gtf_subset_ordr$transcript_id )
trans_seqs_comb <- lapply(exon_seqs_ordr_split, unlist)

library(parallel)
trans_seqs_comb <- mclapply(exon_seqs_ordr_split, unlist, mc.cores = 10)
save(trans_seqs_comb, file = "byCDS/trans_seq.RData")

trans_seq_tri_freq <- lapply(trans_seqs_comb, trinucleotideFrequency)
save(trans_seq_tri_freq, file = "byCDS/trans_seq_tri_freq.RData")

trans_seq_tri_freq_frame1 <- lapply(trans_seqs_comb, function(x) trinucleotideFrequency(x, step = 3))
trans_seq_tri_freq_frame2 <- lapply(trans_seqs_comb, function(x) trinucleotideFrequency(x[-1], step = 3))
trans_seq_tri_freq_frame3 <- lapply(trans_seqs_comb, function(x) trinucleotideFrequency(x[-c(1:2)], step = 3))
save(trans_seq_tri_freq_frame1,trans_seq_tri_freq_frame2,trans_seq_tri_freq_frame3 , file = "byCDS/trans_seq_tri_freq_byFrame.RData")

