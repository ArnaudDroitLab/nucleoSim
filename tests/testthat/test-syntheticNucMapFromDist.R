### Unit tests for syntheticNucMapFromDist.R functions

library(nucleoSim)


test_that("syntheticNucMapFromDist() must generate an error when characters used for as.ratio parameter", {

    message <- "as.ratio must be a logical value"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = 3, fuz.num = 4,
                                            fuz.var = 40,
                                            lin.len = 4, as.ratio = "ICI",
                                            rnd.seed = 15, distr = "Student"),
                        message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when characters used for distr parameter", {

    message <- "\'arg\' should be one of \"Uniform\", \"Normal\", \"Student\""
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = 3, fuz.num = 4,
                                            fuz.var = 40,
                                            lin.len = 4, as.ratio = TRUE,
                                            rnd.seed = 15, distr = "TOTO"),
                        message = message)

})

test_that("syntheticNucMapFromDist() must generate expected results with Normal distribution", {

    obs <- syntheticNucMapFromDist(wp.num = 4, wp.del = 2, wp.var = 3,
                                    fuz.num = 1, fuz.var = 40, lin.len = 4,
                                    as.ratio = FALSE,
                                    rnd.seed = 25, distr = "Normal")
    exp.wp.starts <- c(152, 303)
    exp.wp.nreads <- c(70, 16)
    exp.fuz.starts <- c(390)
    exp.fuz.nreads <- c(57)
    exp.nuc.len <- c(147)

    expect_equal(exp.wp.starts, obs$wp.starts)
    expect_equal(exp.wp.nreads, obs$wp.nreads)
    expect_equal(exp.fuz.starts, obs$fuz.starts)
    expect_equal(exp.fuz.nreads, obs$fuz.nreads)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(86, length(obs$wp.reads))
    expect_equal(57, length(obs$fuz.reads))
    expect_equal(143, length(obs$syn.reads))
    expect_true(is.null(obs$syn.ratio))
    expect_true(is.null(obs$ctr.reads))

})

test_that("syntheticNucMapFromDist() must generate expected results with Student distribution test 01", {

    obs <- syntheticNucMapFromDist(wp.num = 4, wp.del = 1, wp.var = 2,
                                    fuz.num = 1, fuz.var = 40, lin.len = 4,
                                    as.ratio = FALSE,
                                    rnd.seed = 26, distr = "Student")
    exp.wp.starts <- c(1, 152, 303)
    exp.wp.nreads <- c(3, 30, 88)
    exp.fuz.starts <- c(334)
    exp.fuz.nreads <- c(74)
    exp.nuc.len <- 147

    expect_equal(exp.wp.starts, obs$wp.starts)
    expect_equal(exp.wp.nreads, obs$wp.nreads)
    expect_equal(exp.fuz.starts, obs$fuz.starts)
    expect_equal(exp.fuz.nreads, obs$fuz.nreads)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(121, length(obs$wp.reads))
    expect_equal(74, length(obs$fuz.reads))
    expect_equal(195, length(obs$syn.reads))
    expect_true(is.null(obs$syn.ratio))
    expect_true(is.null(obs$ctr.reads))

})

test_that("syntheticNucMapFromDist() must generate expected results with Uniform distribution test 01", {

    obs <- syntheticNucMapFromDist(wp.num = 3, wp.del = 1, nuc.len = 140,
                                    wp.var = 2, fuz.num = 2, fuz.var = 40,
                                    lin.len = 4, as.ratio = FALSE,
                                    rnd.seed = 26, distr = "Uniform")
    exp.wp.starts <- c(1, 145)
    exp.wp.nreads <- c(3, 30)
    exp.fuz.starts <- c(299, 107)
    exp.fuz.nreads <- c(42, 9)
    exp.nuc.len <- 140

    expect_equal(exp.wp.starts, obs$wp.starts)
    expect_equal(exp.wp.nreads, obs$wp.nreads)
    expect_equal(exp.fuz.starts, obs$fuz.starts)
    expect_equal(exp.fuz.nreads, obs$fuz.nreads)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(33, length(obs$wp.reads))
    expect_equal(51, length(obs$fuz.reads))
    expect_equal(84, length(obs$syn.reads))
    expect_true(is.null(obs$syn.ratio))
    expect_true(is.null(obs$ctr.reads))

})

test_that("syntheticNucMapFromDist() must generate expected results with Uniform distribution test 02", {

    obs <- syntheticNucMapFromDist(wp.num = 5, wp.del = 2, wp.var = 3,
                                    fuz.num = 2, fuz.var = 40, lin.len = 4,
                                    as.ratio = TRUE,
                                    rnd.seed = 26, distr = "Uniform")
    exp.wp.starts <- c(1, 152, 454)
    exp.wp.nreads <- c(3, 30, 80)
    exp.fuz.starts <- c(96, 326)
    exp.fuz.nreads <- c(91, 13)
    exp.nuc.len <- 147

    expect_equal(exp.wp.starts, obs$wp.starts)
    expect_equal(exp.wp.nreads, obs$wp.nreads)
    expect_equal(exp.fuz.starts, obs$fuz.starts)
    expect_equal(exp.fuz.nreads, obs$fuz.nreads)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(113, length(obs$wp.reads))
    expect_equal(104, length(obs$fuz.reads))
    expect_equal(217, length(obs$syn.reads))
    expect_equal(664, length(obs$syn.ratio))
    expect_equal(217, length(obs$ctr.reads))

})

test_that("syntheticNucMapFromDist() must generate expected results with Student distribution test 02", {

    obs <- syntheticNucMapFromDist(wp.num = 5, wp.del = 0, max.cover = 15,
                                    wp.var = 23, fuz.num = 2, fuz.var = 43,
                                    lin.len = 12,
                                    as.ratio = FALSE, rnd.seed = 288,
                                    distr = "Student")
    exp.wp.starts <- c(1, 160, 319, 478, 637)
    exp.wp.nreads <- c(3, 11, 4, 11, 3)
    exp.fuz.starts <- c(660, 452)
    exp.fuz.nreads <- c(13, 8)
    exp.nuc.len <- 147

    expect_equal(exp.wp.starts, obs$wp.starts)
    expect_equal(exp.wp.nreads, obs$wp.nreads)
    expect_equal(exp.fuz.starts, obs$fuz.starts)
    expect_equal(exp.fuz.nreads, obs$fuz.nreads)
    expect_equal(exp.nuc.len, obs$nuc.len)
    expect_equal(32, length(obs$wp.reads))
    expect_equal(21, length(obs$fuz.reads))
    expect_equal(53, length(obs$syn.reads))
    expect_equal(0, length(obs$syn.ratio))
    expect_equal(0, length(obs$ctr.reads))

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.num parameter is negative integer", {

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 0,
                                            wp.var = 30, fuz.num = -1,
                                            fuz.var = 40, as.ratio = TRUE,
                                            rnd.seed = 15, distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.var parameter is negative integer", {

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = 1, fuz.num = 4,
                                            fuz.var = -40, as.ratio = TRUE,
                                            rnd.seed = 15, distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when len.var parameter is negative integer", {

    message <- "len.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = 1, fuz.num = 4,
                                            fuz.var = 40, len.var = -1,
                                            as.ratio = TRUE, rnd.seed = 15,
                                            distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when lin.len parameter is negative integer", {

    message <- "len.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = 1, fuz.num = 4,
                                            fuz.var = 40, lin.len = -1,
                                            as.ratio = TRUE, rnd.seed = 15,
                                         distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when max.cover parameter is negative integer",
{

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10,
                                         wp.del = 1, wp.var = 1,
                                         fuz.num = 4, fuz.var = 40,
                                         max.cover = -1,
                                         as.ratio = TRUE, rnd.seed = 15,
                                         distr = "Normal"),
                 message = message)

})


test_that("syntheticNucMapFromDist() must generate an error when wp.del parameter is negative integer",
{

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = -1,
                                            wp.var = 30, fuz.num = 4,
                                            fuz.var = 40, as.ratio = TRUE,
                                            rnd.seed = 15, distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.num parameter is negative integer", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = -1, wp.del = 0,
                                         wp.var = 30, fuz.num = 10, fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.var parameter is negative integer", {

    message <- "wp.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                            wp.var = -1, fuz.num = 4,
                                            fuz.var = 40, as.ratio = TRUE,
                                            rnd.seed = 15, distr = "Normal"),
                    message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.num parameter is character string", {

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = 0,
                                         wp.var = 30, fuz.num = "allo", fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.var parameter is character string", {

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                                         wp.var = 2, fuz.num = 4, fuz.var = "allo", as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when lin.len parameter is character string", {

    message <- "lin.len must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                                         wp.var = 2, fuz.num = 4, fuz.var = 2, lin.len = "allo",
                                         as.ratio = TRUE, rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when max.cover parameter is character string", {

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                                         wp.var = 2, fuz.num = 4, fuz.var = 2, max.cover = "allo",
                                         as.ratio = TRUE, rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.del parameter is character string", {

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = "allo",
                                         wp.var = 30, fuz.num = 4, fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.num parameter is character string", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = "allo",
                                         wp.del = 0, wp.var = 30, fuz.num = 10, fuz.var = 40,
                                         as.ratio = TRUE, rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.var parameter is character string", {

    message <- "wp.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                                         wp.var = "allo", fuz.num = 4, fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.num parameter is a vector of numerics", {

    message <- "fuz.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 0,
                                         wp.var = 30, fuz.num = c(1, 3), fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when fuz.var parameter is a vector of numerics", {

    message <- "fuz.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                         wp.var = 3, fuz.num = 4, fuz.var = c(2, 3), as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when len.var parameter is a vector of numerics", {

    message <- "len.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                         wp.var = 3, fuz.num = 4, fuz.var = 40,
                                         len.var = c(1, 4), as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when lin.len parameter is a vector of numerics", {

    message <- "lin.len must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                         wp.var = 3, fuz.num = 4, fuz.var = 40,
                                         lin.len = c(1, 4), as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when max.cover parameter is a vector of numerics", {

    message <- "max.cover must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                         wp.var = 3, fuz.num = 4, fuz.var = 2,
                                         max.cover = c(1, 2),
                                         as.ratio = TRUE, rnd.seed = 15,
                                         distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.del parameter is a vector of numerics", {

    message <- "wp.del must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = c(1, 3),
                                         wp.var = 30, fuz.num = 4,
                                         fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.num parameter is a vector of numerics", {

    message <- "wp.num must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = c(1, 2),
                                         wp.del = 0, wp.var = 30, fuz.num = 3,
                                         fuz.var = 40, as.ratio = TRUE,
                                         rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when wp.var parameter is a vector of numerics", {

    message <- "wp.var must be a non-negative integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                        wp.var = c(1, 2), fuz.num = 4,
                                        fuz.var = 40, as.ratio = TRUE,
                                        rnd.seed = 15, distr = "Normal"),
                 message = message)

})

test_that("syntheticNucMapFromDist() must generate an error when max.cover parameter is zero", {

    message <- "max.cover must be a positive integer"
    expect_error(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                         wp.var = 1, fuz.num = 4, fuz.var = 40,
                                         max.cover = 0,
                                         as.ratio = TRUE, rnd.seed = 15, distr = "Normal"),
                 message = message)

})
