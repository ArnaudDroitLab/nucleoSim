### Unit tests for syntheticNucReadsFromDist.R functions

library(nucleoSim)

test_that("syntheticNucReadsFromDist() must generate an error when distr parameter is not one of the choice",
{

    message <- "\'arg\' should be one of \"Uniform\", \"Normal\", \"Student\""
    expect_error(syntheticNucReadsFromDist(wp.num = 10, wp.del = 1,
                                           offset = 10, wp.var = 3, fuz.num = 4,
                                           fuz.var = 40, lin.len = 4,
                                           rnd.seed = 15, distr = "TOTO"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must return expected result test 01",
{
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 2,
                                     wp.var = 3, fuz.num = 1, fuz.var = 40,
                                     lin.len = 4, rnd.seed = 125,
                                     offset = 100, distr = "Normal")
    exp.wp <- data.frame(nucleopos = c(175, 326), nreads = c(83,
                                                             13))
    exp.fuz <- data.frame(nucleopos = c(236), nreads = c(61))
    exp.nuc.len <- 147
    exp.dataIP.colnames <- c("chr", "start", "end", "strand",
                             "ID")
    exp.dataIP.chr <- rep("chr_SYNTHETIC", 314)
    exp.paired.colnames <- c("chr", "start", "end", "ID")
    exp.paired.chr <- rep("chr_SYNTHETIC", 157)

    expect_equal(exp.wp, obs$wp)
    expect_equal(exp.fuz, obs$fuz)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(5, length(obs$dataIP))
    expect_equal(exp.dataIP.colnames, colnames(obs$dataIP))
    expect_equal(314, nrow(obs$dataIP))
    expect_equal(exp.dataIP.chr, as.character(obs$dataIP$chr),
                 info = message)
    expect_equal(4, length(obs$paired))
    expect_equal(exp.paired.colnames, colnames(obs$paired))
    expect_equal(157, nrow(obs$paired))
    expect_equal(exp.paired.chr, as.character(obs$paired$chr))

})

test_that("syntheticNucReadsFromDist() must return expected result test 02",
{
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 3,
                                     nuc.len = 144, wp.var = 6, fuz.num = 2,
                                     fuz.var = 60,
                                     lin.len = 4, rnd.seed = 129,
                                     offset = 1000, distr = "Student")
    exp.wp <- data.frame(nucleopos = c(1073), nreads = c(14))
    exp.fuz <- data.frame(nucleopos = c(1416, 1202), nreads = c(23, 87))
    exp.nuc.len <- 144
    exp.dataIP.colnames <- c("chr", "start", "end", "strand", "ID")
    exp.paired.colnames <- c("chr", "start", "end", "ID")
    exp.paired.chr <- rep("chr_SYNTHETIC", 124)

    expect_equal(exp.wp, obs$wp)
    expect_equal(exp.fuz, obs$fuz)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(5, length(obs$dataIP))
    expect_equal(exp.dataIP.colnames, colnames(obs$dataIP))
    expect_equal(248, nrow(obs$dataIP))
    expect_equal(4, length(obs$paired))
    expect_equal(exp.paired.colnames, colnames(obs$paired))
    expect_equal(124, nrow(obs$paired))
    expect_equal(exp.paired.chr, as.character(obs$paired$chr))

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.num parameter is a negative integer",
{

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12,
                                           wp.del = 1, wp.var = 1,
                                           fuz.num = -1, fuz.var = 40,
                                           rnd.seed = 15, offset = 2,
                                           distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.var parameter is a negative integer",
{

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12,
                                           wp.del = 1, wp.var = 1,
                                           fuz.num = 12, fuz.var = -1,
                                           rnd.seed = 15, offset = 2,
                                           distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when lin.len parameter is a negative integer",
{

    message <- "lin.len must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2,
                                           wp.del = 0, wp.var = 0,
                                           fuz.num = 10, fuz.var = 88,
                                           max.cover = 11, lin.len = -1,
                                           rnd.seed = 125, offset = 12,
                                           distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when max.cover parameter is a negative integer",
{

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2,
                                           wp.del = 0, wp.var = 0,
                                           fuz.num = 10, fuz.var = 88,
                                           max.cover = -1, rnd.seed = 125,
                                           offset = 12, distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when nuc.len parameter is a negative integer",
{

              message <- "nuc.len must be a positive integer"
              expect_error(syntheticNucReadsFromDist(wp.num = 2,
                                                     wp.del = 0, wp.var = 0, fuz.num = 10, fuz.var = 88,
                                                     max.cover = 11, nuc.len = -1, rnd.seed = 125, offset = 12,
                                                     distr = "Normal"),
                           message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when offset parameter is a negative integer",
{

    message <- "offset must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                           wp.var = 1, fuz.num = 12, fuz.var = 1, rnd.seed = 15,
                                           offset = -2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when read.len parameter is a negative integer",
{

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12,
                                           wp.del = 1, read.len = -2,
                                           wp.var = 1, fuz.num = 12,
                                           fuz.var = 1, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.del parameter is a negative integer", {

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12, wp.del = -1,
                                           wp.var = 30, fuz.num = 10,
                                           fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.num parameter is a negative integer", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = -1, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.var parameter is a negative integer", {

    message <- "wp.var must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                           wp.var = -1, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.num parameter is a string character", {

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = "10", fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.var parameter is a string character", {

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = "40", rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when lin.len parameter is a string character", {

    message <- "lin.len must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 88, max.cover = 22,
                                           lin.len = "ici", rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when max.cover parameter is a string character",
{

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2,
                                           wp.del = 0, wp.var = 30,
                                           fuz.num = 10, fuz.var = 88,
                                           max.cover = "allo",
                                           rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when nuc.len parameter is a string character", {

    message <- "nuc.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 88, max.cover = 22,
                                           nuc.len = "ici", rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when offset parameter is a string character", {

    message <- "offset must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = "2", distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when read.len parameter is a string character", {

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           read.len = "allo", wp.var = 30, fuz.num = 10, fuz.var = 40,
                                           rnd.seed = 15, offset = 100, distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when rnd.seed parameter is a string character", {

    message <- "rnd.seed must be NULL or an integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = "15",
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.del parameter is a string character", {

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = "0",
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.del parameter is a string character", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = "2", wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.num parameter is a string character", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = "2", wp.del = 0,
                                           wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                           offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.num parameter is a string character", {

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           wp.var = 4, fuz.num = c(1, 3), fuz.var = 40, rnd.seed = 15,
                                           offset = 12, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when fuz.var parameter is a vector of integers", {

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           wp.var = 4, fuz.num = 1, fuz.var = c(1, 3), rnd.seed = 15,
                                           offset = 12, distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when offset parameter is a vector of integers", {

    message <- "offset must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           wp.var = 4, fuz.num = 1, fuz.var = 1, rnd.seed = 15,
                                           offset = c(1, 2), distr = "Student"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when read.len parameter is a vector of integers", {

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           read.len = c(1, 2), wp.var = 4, fuz.num = 1, fuz.var = 1,
                                           rnd.seed = 15, offset = 100, distr = "Student"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when rnd.seed parameter is a vector of integers", {

    message <- "rnd.seed must be NULL or an integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           wp.var = 4, fuz.num = 1,
                                           fuz.var = 2, rnd.seed = c(1, 2),
                                           offset = 12, distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when wp.del parameter is a vector of integers", {

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = c(1, 4),
                                           wp.var = 30, fuz.num = 3,
                                           fuz.var = 40, rnd.seed = 15,
                                           offset = 12, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when wp.num parameter is a vector of integers", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = c(1, 2),
                                           wp.del = 0, wp.var = 30, fuz.num = 3,
                                           fuz.var = 40, rnd.seed = 15,
                                           offset = 3, distr = "Normal"),
                 message = message)
})

test_that("syntheticNucReadsFromDist() must generate an error when wp.var parameter is a vector of integers", {

    message <- "wp.var must be a non-negative integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                           wp.var = c(1, 3), fuz.num = 3,
                                           fuz.var = 40, rnd.seed = 15,
                                           offset = 12, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when max.cover parameter is zero", {

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 0, fuz.num = 10,
                                           fuz.var = 88, max.cover = 0,
                                           rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when nuc.len parameter is zero", {

    message <- "nuc.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                           wp.var = 0, fuz.num = 10,
                                           fuz.var = 88, max.cover = 10,
                                           nuc.len = 0, rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucReadsFromDist() must generate an error when read.len parameter is zero", {

    message <- "read.len must be a positive integer"
    expect_error(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                           read.len = 0, wp.var = 1,
                                           fuz.num = 12, fuz.var = 1,
                                           rnd.seed = 15, offset = 2, distr = "Normal"),
                 message = message)

})
