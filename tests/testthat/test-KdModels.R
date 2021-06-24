data("SampleKdModel", package="scanMiR")

kmers <- data.frame(
  type=c("8mer", "6mer-m8", "7mer-m8", "7mer-a1", "6mer", "6mer-a1",
           "g-bulged 7mer", "wobbled 8mer", "g-bulged 6mer", "wobbled 7mer"),
  log_kd=c(-5083L, -1550L, -4108L, -3749L, -2779L, -1624L, 97L, -1844L, 186L,
           -1216L),
  kmer = c("NNAGCATTAANN", "NNAGCATTCANN", "NNAGCATTACNN", "NNCGCATTAANN",
           "NNCGCATTACNN", "NNAACATTAANN", "NNGCGATTAANN", "NNAGTATTAANN",
           "NNGCGATTACNN", "NNAGTATTACNN")
)

test_that("Kd assignment works", {
  expect_type(assignKdType(kmers$kmer, SampleKdModel)$log_kd, "integer")
  expect_equal(assignKdType(kmers$kmer, SampleKdModel)$log_kd,
               kmers$log_kd, tolerance=200)
})

test_that("Type assignment works", {
  expect_equal(as.character(assignKdType(kmers$kmer, SampleKdModel)$type),
               kmers$type)
})

test_that("KdModel construction works", {
  mod <- getKdModel(kd=dummyKdData(), mirseq="TTAATGCTAATCGTGATAGGGGTT",
                    name="my-miRNA")
  expect_s4_class(mod, "KdModel")
  expect_equal(all(!is.na(mod$mer8) & !is.infinite(mod$mer8)), TRUE)
  expect_equal(all(!is.na(mod$fl) & !is.infinite(mod$fl)), TRUE)
  expect_type(summary(mod), "character")
  expect_gt(mod$cor,0.9)
  expect_gt(cor(mod$mer8, SampleKdModel$mer8),0.9)
  expect_equal(as.character(assignKdType(kmers$kmer, mod)$type),
               kmers$type)
  expect_equal(assignKdType(c("CTAGCATTAAGT","CTAGCATTACGT"), mod)$log_kd,
               c(-5083,-4108), tolerance=250)
})

test_that("KdModelList operations work", {
  expect_s4_class(KdModelList(mod,SampleKdModel), "KdModelList")
  expect_s4_class(ml <- c(mod,SampleKdModel), "KdModelList")
  expect_s4_class(ml[[1]], "KdModel")
  expect_s4_class(ml[2:1], "KdModelList")
  expect_equal(names(ml), c("my-miRNA", "hsa-miR-155-5p"))
})

