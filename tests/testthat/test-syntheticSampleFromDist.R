### Unit tests for syntheticSampleFromDist.R functions

library(nucleoSim)


test_that("test.syntheticNucReadsFromDist_bad_distr", {
    message <- paste0(" syntheticNucReadsFromDist_vector_bad_distr() ",
                      "- A bad value as distr parameter did not generate ",
                      "expected error.")
    expect_error(syntheticNucReadsFromDist(wp.num = 10, wp.del = 1,
                                           offset = 10, wp.var = 3, fuz.num = 4, fuz.var = 40, lin.len = 4,
                                           rnd.seed = 15, distr = "TOTO"), info = message)
})

test_that("test.syntheticNucReadsFromDist_good_result_01", {
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 2,
                                     wp.var = 3, fuz.num = 1, fuz.var = 40, lin.len = 4, rnd.seed = 125,
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
    message <- paste0(" test.syntheticNucReadsFromDist_good_result_01() ",
                      "- syntheticNucReadsFromDist did not generate ", "expected values")
    expect_equal(exp.wp, obs$wp, info = message)
    expect_equal(exp.fuz, obs$fuz, info = message)
    expect_equal(exp.nuc.len, obs$nuc.len, info = message)
    expect_equal(5, length(obs$dataIP), info = message)
    expect_equal(exp.dataIP.colnames, colnames(obs$dataIP), info = message)
    expect_equal(314, nrow(obs$dataIP), info = message)
    expect_equal(exp.dataIP.chr, as.character(obs$dataIP$chr),
                 info = message)
    expect_equal(4, length(obs$paired), info = message)
    expect_equal(exp.paired.colnames, colnames(obs$paired), info = message)
    expect_equal(157, nrow(obs$paired), info = message)
    expect_equal(exp.paired.chr, as.character(obs$paired$chr),
                 info = message)
})

test_that("test.syntheticNucReadsFromDist_good_result_02", {
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 3,
                                     nuc.len = 144, wp.var = 6, fuz.num = 2, fuz.var = 60,
                                     lin.len = 4, rnd.seed = 129, offset = 1000, distr = "Student")
    exp.wp <- data.frame(nucleopos = c(1073), nreads = c(14))
    exp.fuz <- data.frame(nucleopos = c(1416, 1202), nreads = c(23,
                                                                87))
    exp.nuc.len <- 144
    exp.dataIP.colnames <- c("chr", "start", "end", "strand",
                             "ID")
    exp.paired.colnames <- c("chr", "start", "end", "ID")
    exp.paired.chr <- rep("chr_SYNTHETIC", 124)
    message <- paste0(" test.syntheticNucReadsFromDist_good_result_02() ",
                      "- syntheticNucReadsFromDist did not generate ", "expected values")
    expect_equal(exp.wp, obs$wp, info = message)
    expect_equal(exp.fuz, obs$fuz, info = message)
    expect_equal(exp.nuc.len, obs$nuc.len, info = message)
    expect_equal(5, length(obs$dataIP), info = message)
    expect_equal(exp.dataIP.colnames, colnames(obs$dataIP), info = message)
    expect_equal(248, nrow(obs$dataIP), info = message)
    expect_equal(4, length(obs$paired), info = message)
    expect_equal(exp.paired.colnames, colnames(obs$paired), info = message)
    expect_equal(124, nrow(obs$paired), info = message)
    expect_equal(exp.paired.chr, as.character(obs$paired$chr),
                 info = message)
})

test_that("test.syntheticNucReadsFromDist_negative_fuz_num",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12,
                                                        wp.del = 1, wp.var = 1, fuz.num = -1, fuz.var = 40,
                                                        rnd.seed = 15, offset = 2, distr = "Normal"), error = conditionMessage)
              exp <- "fuz.num must be a non-negative integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_fuz_num() ",
                                "- A negative integer as fuz.num parameter did not generate ",
                                "expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_fuz_var",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12,
                                                        wp.del = 1, wp.var = 1, fuz.num = 12, fuz.var = -1,
                                                        rnd.seed = 15, offset = 2, distr = "Normal"), error = conditionMessage)
              exp <- "fuz.var must be a non-negative integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_fuz_var() ",
                                "- A negative integer as fuz.var parameter did not generate ",
                                "expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_lin_len",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2,
                                                        wp.del = 0, wp.var = 0, fuz.num = 10, fuz.var = 88,
                                                        max.cover = 11, lin.len = -1, rnd.seed = 125, offset = 12,
                                                        distr = "Normal"), error = conditionMessage)
              exp <- "lin.len must be a non-negative integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_lin_len() ",
                                "- Negative value as lin.len parameter did not ",
                                "generate expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_max_cover",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2,
                                                        wp.del = 0, wp.var = 0, fuz.num = 10, fuz.var = 88,
                                                        max.cover = -1, rnd.seed = 125, offset = 12, distr = "Normal"),
                              error = conditionMessage)
              exp <- "max.cover must be a positive integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_max_cover() ",
                                "- Negative value as max.cover parameter did not ",
                                "generate expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_nuc_len",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2,
                                                        wp.del = 0, wp.var = 0, fuz.num = 10, fuz.var = 88,
                                                        max.cover = 11, nuc.len = -1, rnd.seed = 125, offset = 12,
                                                        distr = "Normal"), error = conditionMessage)
              exp <- "nuc.len must be a positive integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_nuc_len() ",
                                "- Negative value as nuc.len parameter did not ",
                                "generate expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_offset", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                              wp.var = 1, fuz.num = 12, fuz.var = 1, rnd.seed = 15,
                                              offset = -2, distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_offset() ",
                      "- A negative integer as offset  parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_negative_read_len",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12,
                                                        wp.del = 1, read.len = -2, wp.var = 1, fuz.num = 12,
                                                        fuz.var = 1, rnd.seed = 15, offset = 2, distr = "Normal"),
                              error = conditionMessage)
              exp <- "read.len must be a positive integer"
              message <- paste0(" test.syntheticNucReadsFromDist_negative_read_len() ",
                                "- A negative integer as read.len parameter did not generate ",
                                "expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_negative_wp_del", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = -1,
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_del() ",
                      "- A negative integer as wp.del parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_negative_wp_num", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = -1, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_num() ",
                      "- A negative integer as wp.num parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_negative_wp_var", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                              wp.var = -1, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_var() ",
                      "- A negative integer as wp.var parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_fuz_num", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = "10", fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_fuz_num() ",
                      "- A string as fuz.num parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_fuz_var", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = "40", rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_fuz_var() ",
                      "- A string as fuz.var parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_lin_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 88, max.cover = 22,
                                              lin.len = "ici", rnd.seed = 15, offset = 2, distr = "Normal"),
                    error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_lin_len() ",
                      "- A string as lin.len parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_max_cover",
          {
              obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2,
                                                        wp.del = 0, wp.var = 30, fuz.num = 10, fuz.var = 88,
                                                        max.cover = "allo", rnd.seed = 15, offset = 2, distr = "Normal"),
                              error = conditionMessage)
              exp <- "max.cover must be a positive integer"
              message <- paste0(" test.syntheticNucReadsFromDist_string_max_cover() ",
                                "- A string as max.cover parameter did not generate ",
                                "expected error.")
              expect_equal(exp, obs, info = message)
          })

test_that("test.syntheticNucReadsFromDist_string_nuc_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 88, max.cover = 22,
                                              nuc.len = "ici", rnd.seed = 15, offset = 2, distr = "Normal"),
                    error = conditionMessage)
    exp <- "nuc.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_nuc_len() ",
                      "- A string as nuc.len parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_offset", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = "2", distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_offset() ",
                      "- A string as offset parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_read_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              read.len = "allo", wp.var = 30, fuz.num = 10, fuz.var = 40,
                                              rnd.seed = 15, offset = 100, distr = "Normal"), error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_read_len() ",
                      "- A string as read.len parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_rnd_seed", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = "15",
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_rnd_seed() ",
                      "- A string as rnd.seed parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_wp_del", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = "0",
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_del() ",
                      "- A string as wp.del parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_wp_num", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = "2", wp.del = 0,
                                              wp.var = 30, fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_num() ",
                      "- A string as wp.num parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_string_wp_var", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = "30", fuz.num = 10, fuz.var = 40, rnd.seed = 15,
                                              offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_var() ",
                      "- A string as wp.var parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_fuz_num", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              wp.var = 4, fuz.num = c(1, 3), fuz.var = 40, rnd.seed = 15,
                                              offset = 12, distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_fuz_num() ",
                      "- A vector as fuz.num parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_fuz_var", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              wp.var = 4, fuz.num = 1, fuz.var = c(1, 3), rnd.seed = 15,
                                              offset = 12, distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_fuz_var() ",
                      "- A vector as fuz.var parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_offset", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              wp.var = 4, fuz.num = 1, fuz.var = 1, rnd.seed = 15,
                                              offset = c(1, 2), distr = "Student"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_offset() ",
                      "- A vector as offset parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_read_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              read.len = c(1, 2), wp.var = 4, fuz.num = 1, fuz.var = 1,
                                              rnd.seed = 15, offset = 100, distr = "Student"), error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_read_len() ",
                      "- A vector as _read.len parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_rnd_seed", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              wp.var = 4, fuz.num = 1, fuz.var = 2, rnd.seed = c(1,
                                                                                                 2), offset = 12, distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_rnd_seed() ",
                      "- A vector as rnd.seed parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_wp_del", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = c(1,
                                                                      4), wp.var = 30, fuz.num = 3, fuz.var = 40, rnd.seed = 15,
                                              offset = 12, distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_wp_del() ",
                      "- A vector as wp.del parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_wp_num", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = c(1, 2),
                                              wp.del = 0, wp.var = 30, fuz.num = 3, fuz.var = 40, rnd.seed = 15,
                                              offset = 3, distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucReadsFromDist_vector_wp_num() ",
                      "- A vector as wp.num parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_vector_wp_var", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                              wp.var = c(1, 3), fuz.num = 3, fuz.var = 40, rnd.seed = 15,
                                              offset = 12, distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_wp_var() ",
                      "- A vector as wp.var parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_zero_max_cover", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 0, fuz.num = 10, fuz.var = 88, max.cover = 0,
                                              rnd.seed = 15, offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "max.cover must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_zero_max_cover() ",
                      "- Zero as max.cover parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_zero_nuc_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                              wp.var = 0, fuz.num = 10, fuz.var = 88, max.cover = 10,
                                              nuc.len = 0, rnd.seed = 15, offset = 2, distr = "Normal"),
                    error = conditionMessage)
    exp <- "nuc.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_zero_nuc_len() ",
                      "- Zero as nuc.len parameter did not generate ", "expected error.")
    expect_equal(exp, obs, info = message)
})

test_that("test.syntheticNucReadsFromDist_zero_read_len", {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                              read.len = 0, wp.var = 1, fuz.num = 12, fuz.var = 1,
                                              rnd.seed = 15, offset = 2, distr = "Normal"), error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromDist_zero_read_len() ",
                      "- A negative integer as read.len parameter did not generate ",
                      "expected error.")
    expect_equal(exp, obs, info = message)
})
