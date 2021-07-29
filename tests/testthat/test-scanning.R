data("SampleKdModel", package="scanMiR")
seq1 <- DNAStringSet(c(seq = paste(c(rep("N", 20), "AGCATTAA", rep("N", 20)),
                                   collapse = "")))

gr <- findSeedMatches(seq1, SampleKdModel, verbose = FALSE, ret = "GRanges")
df <- findSeedMatches(seq1, SampleKdModel, verbose = FALSE, ret = "data.frame")
agg <- findSeedMatches(seq1, SampleKdModel, verbose = FALSE, ret = "aggregated")

test_that("Scan returns the right format", {
  expect_s4_class(gr, "GRanges")
  expect_s3_class(df, "data.frame")
  expect_s3_class(agg, "data.frame")
})

test_that("Scan returns the right match", {
  expect_identical(as.character(seqnames(gr)), "seq")
  expect_identical(ranges(gr), IRanges(21, 28))
  expect_identical(as.character(gr$type), "8mer")
})

seq2 <- DNAStringSet(c(seq = paste(c(rep("N", 5), "AGCATTAA", rep("N", 7),
                                     "AGCATTAA", rep("N", 20)), collapse = "")))
gr_s0 <- findSeedMatches(seq2, SampleKdModel, verbose = FALSE, shadow = 0L)
gr_s15 <- findSeedMatches(seq2, SampleKdModel, verbose = FALSE, shadow = 15L)

test_that("Shadow works", {
  expect_identical(length(gr_s0), 2L)
  expect_identical(length(gr_s15), 1L)
})

seq_nc <- DNAStringSet(c(seq = paste(c(rep("N", 5), "AGTATTAA", rep("N", 7),
                                       "AGCATTAA", rep("N", 20)),
                                     collapse = "")))
gr_withNC <- findSeedMatches(seq_nc, SampleKdModel, verbose = FALSE,
                             onlyCanonical = FALSE)
gr_withoutNC <- findSeedMatches(seq_nc, SampleKdModel, verbose = FALSE,
                                onlyCanonical = TRUE)

test_that("OnlyCanonical works", {
  expect_identical(length(gr_withNC), 2L)
  expect_identical(length(gr_withoutNC), 1L)
})

gr_ms <- findSeedMatches(seq1, SampleKdModel, verbose = FALSE, ret = "GRanges",
                         keepMatchSeq = TRUE)
test_that("KeepMatchSeq works", {
  expect_identical(as.character(gr_ms$sequence$seq), "NNAGCATTAANN")
})

gr_md_0 <- findSeedMatches(seq_nc, SampleKdModel, verbose = FALSE,
                           onlyCanonical = FALSE, minDist = 0L)
gr_md_10 <- findSeedMatches(seq_nc, SampleKdModel, verbose = FALSE,
                            onlyCanonical = FALSE, minDist = 10L)
test_that("minDist works", {
  expect_identical(length(gr_md_0), 2L)
  expect_identical(length(gr_md_10), 1L)
  expect_identical(as.character(gr_md_10$type), "8mer")
})

char_seq <- setNames(as.character(seq1$seq), "seq")
gr_char <- findSeedMatches(char_seq, SampleKdModel, verbose = FALSE,
                           ret = "GRanges")
test_that("Scan works for character sequence", {
  expect_identical(gr_char, gr)
})

char_seed <- SampleKdModel$canonical.seed
gr_char_seed <- findSeedMatches(seq1, char_seed, verbose = FALSE,
                                ret = "GRanges")
test_that("Scan works for character seed", {
  expect_identical(ranges(gr), ranges(gr_char_seed))
  expect_identical(gr$type, gr_char_seed$type)
})

test_that("Aggregation works", {
  expect_identical(nrow(agg), 1L)
  expect_identical(names(agg)[2], "repression")
})

seq_ORF <- seq2
mcols(seq_ORF)$ORF.length <- 15
gr_ORF <- findSeedMatches(seq_ORF, SampleKdModel, verbose = FALSE)
df_ORF <- findSeedMatches(seq_ORF, SampleKdModel, verbose = FALSE,
                          ret = "data.frame")
agg_ORF <- findSeedMatches(seq_ORF, SampleKdModel, verbose = FALSE,
                          ret = "aggregated")
test_that("Scanning ORF region works", {
  expect_identical(as.logical(gr_ORF$ORF), c(TRUE, FALSE))
  expect_s3_class(metadata(gr_ORF)$tx_info, "data.frame")
  expect_identical(names(metadata(gr_ORF)$tx_info), c("ORF.length", "length",
                                                      "UTR.length"))
  expect_identical(as.numeric(metadata(gr_ORF)$tx_info), c(15,48,33))
  expect_identical(df_ORF$ORF, c(TRUE, FALSE))
  expect_identical(attr(df_ORF, "tx_info"), metadata(gr_ORF)$tx_info)
  expect_identical(attr(agg_ORF, "tx_info"), metadata(gr_ORF)$tx_info)
  expect_identical(agg_ORF$ORF.canonical, 1L)
})

x <- c("AACACTCCAG","GACACTCCGC","GTACTCCAT","ACGTACGTAC")
matchTypes <- as.character(getMatchTypes(x, seed="ACACTCCA"))
test_that("Correct match types are returned", {
  expect_identical(matchTypes, c("8mer", "7mer-m8", "6mer-a1", "non-canonical"))
})

gr_OR_0 <- removeOverlappingRanges(gr_withNC, minDist = 0L)
gr_OR_10 <- removeOverlappingRanges(gr_withNC, minDist = 10L)
test_that("Correctly removes overlapping ranges", {
  expect_identical(length(gr_OR_0), 2L)
  expect_identical(length(gr_OR_10), 1L)
  expect_identical(as.character(gr_OR_10$type), "8mer")
})

test_that("3p alignment works correctly", {
  al3p <- get3pAlignment(seqs="TGTTTAAGAACAACAAAATCACCAT",
                         mirseq="TGGAAGACTAGTGATTTTGTTGTT", siteType="8mer")
  expect_equal(as.numeric(al3p[1,1:4]), c(2,3,0,14))
  expect_identical(as.character(al3p$note), "TDMD?")
  al3p <- get3pAlignment( "GCTAAGTTCTCCCAACAACATGAA", "TAGGTAGTTTCATGTTGTTGGG",
                          siteType="8mer")
  expect_equal(as.numeric(al3p[1,1:4]), c(0,0,0,14))
  expect_identical(as.character(al3p$note), "Slicing")
})

test_that("viewTargetAlignment works", {
  expect_s3_class(viewTargetAlignment(gr_OR_0[1], SampleKdModel, seq_nc,
                                      outputType="ggplot"), "ggplot")
  al <- viewTargetAlignment(gr_OR_0[1], SampleKdModel, seq_nc,
                            outputType="data.frame")
  expect_identical(al$alignment[2],
                   "                   |||-|        ||||||||        ")
})

