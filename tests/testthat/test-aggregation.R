matches <- GRanges("seq",
                   IRanges(start = c(1,10), width = 7),
                   type = c("8mer", "7mer"),
                   log_kd = c(-5000, -3000),
                   p3.score = 8
)

agg <- aggregateMatches(matches, p3 = 0, coef_utr = 0, coef_orf = 0)

test_that("Aggregation format is correct", {
  expect_s3_class(agg, "data.frame")
  expect_identical(names(agg), c("transcript", "repression", "8mer", "7mer", 
                                 "6mer", "non-canonical"))
})

test_that("Predicted repression is correct", {
  expect_equal(agg$repression, -1.089705, tolerance = 1e-6)
})

test_that("Site info is correct", {
  expect_identical(agg$`8mer`, 1L)
  expect_identical(agg$`7mer`, 1L)
  expect_identical(agg$`6mer`, 0L)
  expect_identical(agg$`non-canonical`, 0L)
})

matches_df <- as.data.frame(matches)
names(matches_df)[1] <- "transcript"
agg_from_df <- aggregateMatches(matches_df, p3 = 0, coef_utr = 0, coef_orf = 0)

test_that("Aggregation also works for data.frame", {
  expect_identical(agg_from_df, agg)
})

matches_ORF <- matches
matches_ORF$ORF <- c(TRUE, FALSE)
agg_ORF <- aggregateMatches(matches_ORF, p3 = 0, coef_utr = 0, coef_orf = 0)
test_that("Predicted repression with ORF site is correct", {
  expect_equal(agg_ORF$repression, -0.6030502, tolerance = 1e-6)
  expect_identical(agg_ORF$`8mer`, 0L)
  expect_identical(agg_ORF$`7mer`, 1L)
  expect_identical(agg_ORF$`6mer`, 0L)
  expect_identical(agg_ORF$`non-canonical`, 0L)  
  expect_identical(agg_ORF$`ORF.canonical`, 1L)
})

agg_3p <- aggregateMatches(matches, p3 = 1, coef_utr = 0, coef_orf = 0)

test_that("Predicted repression with 3p score is correct", {
  expect_equal(agg_3p$repression, -2.14542, tolerance = 1e-6)
})

matches_noKD <- matches
matches_noKD$log_kd <- NULL
agg_noKD <- aggregateMatches(matches_noKD)

test_that("Site info is returned if no KDs are given", {
  expect_identical(agg_noKD, agg[,-2])
})

matches_miRNA <- c(matches_ORF, matches_ORF)
matches_miRNA$miRNA <- factor(rep(c("miRNA1", "miRNA2"), each = 2))
agg_miRNA <- aggregateMatches(matches_miRNA, p3 = 0, coef_utr = 0, coef_orf = 0)
test_that("Aggregation works with several miRNAs", {
  expect_identical(agg_miRNA$miRNA, factor(c("miRNA1", "miRNA2")))
  expect_identical(agg_miRNA[,-1], rbind(agg_ORF, agg_ORF))
})

matches_coef <- GRanges(seqnames = c("seq1", "seq2"), 
                        ranges(matches),
                        type = matches$type,
                        log_kd = matches$log_kd,
                        p3.score = matches$p3.score,
                        ORF = matches_ORF$ORF,
                        miRNA = "miRNA1"
                        )
metadata(matches_coef)$tx_info <- data.frame(ORF.length = c(10, 20),
                                              UTR.length = c(10, 20),
                                              row.names = c("seq1", "seq2"))
agg_coef <- aggregateMatches(matches_coef, p3 = 0, coef_utr = 1, coef_orf = 1)
test_that("Predicted repression with coef_utr and coef_orf is correct", {
  expect_equal(agg_coef$repression, c(-0.338184, -0.898957), tolerance = 1e-6)
})