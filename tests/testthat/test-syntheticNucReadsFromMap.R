### Unit tests for syntheticNucReadsFromMap.R functions

library(nucleoSim)


NUCLEO_MAP <- syntheticNucMapFromDist(wp.num=3, wp.del=1, wp.var=12,
                                      fuz.num=1, fuz.var=44, max.cover=65,
                                      nuc.len=147, len.var=2, lin.len=40,
                                      rnd.seed=155, distr="Uniform")


test_that("syntheticNucReadsFromMap( must return the expected results", {

    obs <- syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                    read.len = 40, offset = 1000)
    exp.wp <- data.frame(nucleopos = c(1262, 1449), nreads = c(29, 50))
    exp.fuz <- data.frame(nucleopos = c(1481), nreads = c(25))
    exp.nuc.len <- 147
    exp.dataIP.colnames <- c("chr", "start", "end", "strand", "ID")
    exp.dataIP.chr <- rep("chr_SYNTHETIC", 208)
    exp.paired.colnames <- c("chr", "start", "end", "ID")
    exp.paired.chr <- rep("chr_SYNTHETIC", 104)

    expect_equal(exp.wp, obs$wp)
    expect_equal(exp.fuz, obs$fuz)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(5, length(obs$dataIP))
    expect_equal(exp.dataIP.colnames, colnames(obs$dataIP))
    expect_equal(208, nrow(obs$dataIP))
    expect_equal(exp.dataIP.chr, as.character(obs$dataIP$chr))
    expect_equal(4, length(obs$paired))
    expect_equal(exp.paired.colnames, colnames(obs$paired))
    expect_equal(104, nrow(obs$paired))
    expect_equal(exp.paired.chr, as.character(obs$paired$chr))

})

test_that("syntheticNucReadsFromMap() must generate an error when syntheticNucMap parameter is an integer",
{

    message <- "syntheticNucMap must be an object of class 'syntheticNucMap'"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = 22,
                                          read.len = 40, offset = 12),
                 message = message)

})

test_that("syntheticNucReadsFromMap() must generate an error when offset parameter is negative integer",
{

    message <- "offset must be a non-negative integer"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                          read.len = 20, offset = -1),
                 message = message)

})

test_that("syntheticNucReadsFromMap() must generate an error when read.len parameter is a negative integer",
{

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                          read.len = -1, offset = 12),
                 message = message)

})

test_that("syntheticNucReadsFromMap() must generate an error when offset parameter is character string", {

    message <- "offset must be a non-negative integer"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                          read.len = 30, offset = "testing"),
                 message = message)

})

test_that("syntheticNucReadsFromMap() must generate an error when read.len parameter is character string",
{

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                          read.len = "igloo", offset = 12),
                 message = message)
})

test_that("syntheticNucReadsFromMap() must generate an error when syntheticNucMap parameter is character string",
{

    message <- "syntheticNucMap must be an object of class 'syntheticNucMap'"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = "canada",
                                          read.len = 40, offset = 12),
                 message = message)
})

test_that("syntheticNucReadsFromMap() must generate an error when read.len parameter is zero",
{

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                          read.len = 0, offset = 12),
                 message = message)

})
