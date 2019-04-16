###################################################
# Created by Astrid Deschenes
# 2015-07-23
###################################################

###################################################
## Test the syntheticNucMapFromDist function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "nucleoSim" )
}

### }}}


###########################################################
## Test the syntheticNucMapFromDist() function parameters
###########################################################

################################
## wp.num
################################

## Test the result when a string is passed as wp.num parameter
test.syntheticNucMapFromDist_string_wp_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = "allo", wp.del = 0,
                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                            as.ratio = TRUE, rnd.seed = 15,
                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_wp_num() ",
                      "- A string as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.num parameter
test.syntheticNucMapFromDist_negative_wp_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = -1, wp.del = 0,
                        wp.var = 30, fuz.num = 10, fuz.var = 40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_wp_num() ",
                "- A negative integer as wp.num parameter did not generated ",
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.num parameter
test.syntheticNucMapFromDist_vector_wp_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = c(1, 2), wp.del = 0,
                        wp.var = 30, fuz.num = 3, fuz.var = 40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_wp_num() ",
                        "- A vector as wp.num parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## fuz.num
################################

## Test the result when a string is passed as fuz.num parameter
test.syntheticNucMapFromDist_string_fuz_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 0,
                        wp.var = 30, fuz.num = "allo", fuz.var = 40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_fuz_num() ",
                        "- A string as fuz.num parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.num parameter
test.syntheticNucMapFromDist_negative_fuz_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 0,
                            wp.var = 30, fuz.num = -1, fuz.var = 40,
                            as.ratio = TRUE, rnd.seed = 15,
                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_fuz_num() ",
                "- A negative integer as fuz.num parameter did not generated ",
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.num parameter
test.syntheticNucMapFromDist_vector_fuz_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 0,
                wp.var = 30, fuz.num = c(1, 3), fuz.var = 40,
                as.ratio = TRUE, rnd.seed = 15,
                distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_fuz_num() ",
                "- A vector as fuz.num parameter did not generated ",
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## wp.del
################################

## Test the result when a string is passed as wp.del parameter
test.syntheticNucMapFromDist_string_wp_del <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = "allo",
                    wp.var = 30, fuz.num = 4, fuz.var = 40,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_wp_del() ",
                        "- A string as wp.del parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.del parameter
test.syntheticNucMapFromDist_negative_wp_del <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = -1,
                        wp.var = 30, fuz.num = 4, fuz.var = 40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_wp_del() ",
            "- A negative integer as wp.del parameter did not generated ",
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.del parameter
test.syntheticNucMapFromDist_vector_wp_del<- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = c(1, 3),
                    wp.var = 30, fuz.num = 4, fuz.var = 40,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_wp_del() ",
                        "- A vector as wp.del parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## wp.var
################################

## Test the result when a string is passed as wp.var parameter
test.syntheticNucMapFromDist_string_wp_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                            wp.var = "allo", fuz.num = 4, fuz.var = 40,
                            as.ratio = TRUE, rnd.seed = 15,
                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_wp_var() ",
                        "- A string as wp.var parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.var parameter
test.syntheticNucMapFromDist_negative_wp_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                        wp.var = -1, fuz.num = 4, fuz.var = 40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_wp_var() ",
                "- A negative integer as wp.var parameter did not generated ",
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.var parameter
test.syntheticNucMapFromDist_vector_wp_var<- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = c(1, 2), fuz.num = 4, fuz.var = 40,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_wp_var() ",
                        "- A vector as wp.var parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## fuz.var
################################

## Test the result when a string is passed as fuz.var parameter
test.syntheticNucMapFromDist_string_fuz_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                    wp.var = 2, fuz.num = 4, fuz.var = "allo",
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_fuz_var() ",
                        "- A string as fuz.var parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.var parameter
test.syntheticNucMapFromDist_negative_fuz_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                        wp.var = 1, fuz.num = 4, fuz.var = -40,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_fuz_var() ",
        "- A negative integer as fuz.var parameter did not generated ",
        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.var parameter
test.syntheticNucMapFromDist_vector_fuz_var<- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                            wp.var = 3, fuz.num = 4, fuz.var = c(2, 3),
                            as.ratio = TRUE, rnd.seed = 15,
                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_fuz_var() ",
                        "- A vector as fuz.var parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## max.cover
################################

## Test the result when a string is passed as max.cover parameter
test.syntheticNucMapFromDist_string_max_cover <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                wp.var = 2, fuz.num = 4, fuz.var = 2, max.cover = "allo",
                as.ratio = TRUE, rnd.seed = 15,
                distr = "Normal"), error = conditionMessage)
    exp <- "max.cover must be a positive integer"
    message <- paste0(" syntheticNucMapFromDist_string_max_cover() ",
                      "- A string as max.cover parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a zero is passed as max.cover parameter
test.syntheticNucMapFromDist_zero_max_cover <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = 1, fuz.num = 4, fuz.var = 40, max.cover = 0,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "max.cover must be a positive integer"
    message <- paste0(" syntheticNucMapFromDist_zero_max_cover() ",
            "- A zero integer as max.cover parameter did not generated ",
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as max.cover parameter
test.syntheticNucMapFromDist_negative_max_cover <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = 1, fuz.num = 4, fuz.var = 40, max.cover = -1,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "max.cover must be a positive integer"
    message <- paste0(" syntheticNucMapFromDist_negative_max_cover() ",
            "- A negative integer as max.cover parameter did not generated ",
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as max.cover parameter
test.syntheticNucMapFromDist_vector_max_cover <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
            wp.var = 3, fuz.num = 4, fuz.var = 2, max.cover = c(1, 2),
            as.ratio = TRUE,  rnd.seed = 15,
            distr = "Normal"), error = conditionMessage)
    exp <- "max.cover must be a positive integer"
    message <- paste0(" syntheticNucMapFromDist_vector_max_cover() ",
                        "- A vector as max.cover parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a string is passed as lin.len parameter
test.syntheticNucMapFromDist_string_lin_len <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                        wp.var = 2, fuz.num = 4, fuz.var = 2, lin.len = "allo",
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_lin_len() ",
                        "- A string as lin.len parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## lin.len
################################

## Test the result when a negative is passed as lin.len parameter
test.syntheticNucMapFromDist_negative_lin_len <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                        wp.var = 1, fuz.num = 4, fuz.var = 40, lin.len = -1,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_negative_lin_len() ",
            "- A negative integer as lin.len parameter did not generated ",
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as lin.len parameter
test.syntheticNucMapFromDist_vector_lin_len <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = 3, fuz.num = 4, fuz.var = 40, lin.len = c(1, 4),
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_lin_len() ",
                        "- A vector as lin.len parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## len.var
################################

## Test the result when a negative is passed as len.var parameter
test.syntheticNucMapFromDist_negative_len_var<- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = 1, fuz.num = 4, fuz.var = 40, len.var = -1,
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "len.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucMapFromDist_negative_len_var() ",
                        "- A negative integer as len.var parameter did not ",
                        "generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as len.var parameter
test.syntheticNucMapFromDist_vector_len_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                    wp.var = 3, fuz.num = 4, fuz.var = 40, len.var = c(1, 4),
                    as.ratio = TRUE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "len.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucMapFromDist_vector_len_var() ",
                        "- A vector as len.var parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## as.ratio
################################

## Test the result when bad value is passed as distr parameter
test.syntheticNucMapFromDist_bad_as_ratio <- function() {
    message <- paste0(" syntheticNucMapFromDist_vector_bad_as_ratio() ",
                      "- A bad value as as.ratio parameter did not generated ",
                      "expected error.")

    checkException(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                                           wp.var = 3, fuz.num = 4,
                                           fuz.var = 40, lin.len = 4,
                                           as.ratio = "ICI", rnd.seed = 15,
                                           distr = "Student"), msg = message)
}


################################
## distr
################################

## Test the result when bad value is passed as distr parameter
test.syntheticNucMapFromDist_bad_distr <- function() {
    message <- paste0(" syntheticNucMapFromDist_vector_bad_distr() ",
                        "- A bad value as distr parameter did not generated ",
                        "expected error.")

    checkException(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                        wp.var = 3, fuz.num = 4, fuz.var = 40, lin.len = 4,
                        as.ratio = TRUE, rnd.seed = 15,
                        distr = "TOTO"), msg = message)
}


################################
## Good results
################################


## Test the result when as.ratio is FALSE
test.syntheticNucMapFromDist_good_result_01 <- function() {
    obs <- syntheticNucMapFromDist(wp.num = 4, wp.del = 2,
                wp.var = 3, fuz.num = 1, fuz.var = 40, lin.len = 4,
                as.ratio = FALSE, rnd.seed = 25,
                distr = "Normal")

    exp.wp.starts  <- c(152, 303)
    exp.wp.nreads  <- c(70, 16)
    exp.fuz.starts <- c(390)
    exp.fuz.nreads <- c(57)
    exp.nuc.len    <- c(147)

    message     <- paste0(" test.syntheticNucMapFromDist_good_result_01() ",
                          "- syntheticNucMapFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp.starts, exp.wp.starts, msg = message)
    checkEqualsNumeric(obs$wp.nreads, exp.wp.nreads, msg = message)
    checkEqualsNumeric(obs$fuz.starts, exp.fuz.starts, msg = message)
    checkEqualsNumeric(obs$fuz.nreads, exp.fuz.nreads, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$wp.reads), 86, msg = message)
    checkEqualsNumeric(length(obs$fuz.reads), 57, msg = message)
    checkEqualsNumeric(length(obs$syn.reads), 143, msg = message)
    checkTrue(is.null(obs$syn.ratio), msg = message)
    checkTrue(is.null(obs$ctr.reads), msg = message)
}

## Test the result when as.ratio is FALSE and Student distribution
test.syntheticNucMapFromDist_good_result_02 <- function() {
    obs <- syntheticNucMapFromDist(wp.num = 4, wp.del = 1,
                                   wp.var = 2, fuz.num = 1, fuz.var = 40,
                                   lin.len = 4,
                                   as.ratio = FALSE, rnd.seed = 26,
                                   distr = "Student")

    exp.wp.starts  <- c(1, 152, 303)
    exp.wp.nreads  <- c(3, 30, 88)
    exp.fuz.starts <- c(334)
    exp.fuz.nreads <- c(74)
    exp.nuc.len    <- 147

    message     <- paste0(" test.syntheticNucMapFromDist_good_result_02() ",
                          "- syntheticNucMapFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp.starts, exp.wp.starts, msg = message)
    checkEqualsNumeric(obs$wp.nreads, exp.wp.nreads, msg = message)
    checkEqualsNumeric(obs$fuz.starts, exp.fuz.starts, msg = message)
    checkEqualsNumeric(obs$fuz.nreads, exp.fuz.nreads, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$wp.reads), 121, msg = message)
    checkEqualsNumeric(length(obs$fuz.reads), 74, msg = message)
    checkEqualsNumeric(length(obs$syn.reads), 195, msg = message)
    checkTrue(is.null(obs$syn.ratio), msg = message)
    checkTrue(is.null(obs$ctr.reads), msg = message)
}

## Test the result when as.ratio is FALSE and Uniform distribution
test.syntheticNucMapFromDist_good_result_03 <- function() {
    obs <- syntheticNucMapFromDist(wp.num = 3, wp.del = 1, nuc.len = 140,
                                    wp.var = 2, fuz.num = 2, fuz.var = 40,
                                    lin.len = 4, as.ratio = FALSE,
                                    rnd.seed = 26, distr = "Uniform")

    exp.wp.starts  <- c(1, 145)
    exp.wp.nreads  <- c(3, 30)
    exp.fuz.starts <- c(299, 107)
    exp.fuz.nreads <- c(42, 9)
    exp.nuc.len    <- 140

    message     <- paste0(" test.syntheticNucMapFromDist_good_result_03() ",
                          "- syntheticNucMapFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp.starts, exp.wp.starts, msg = message)
    checkEqualsNumeric(obs$wp.nreads, exp.wp.nreads, msg = message)
    checkEqualsNumeric(obs$fuz.starts, exp.fuz.starts, msg = message)
    checkEqualsNumeric(obs$fuz.nreads, exp.fuz.nreads, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$wp.reads), 33, msg = message)
    checkEqualsNumeric(length(obs$fuz.reads), 51, msg = message)
    checkEqualsNumeric(length(obs$syn.reads), 84, msg = message)
    checkTrue(is.null(obs$syn.ratio), msg = message)
    checkTrue(is.null(obs$ctr.reads), msg = message)
}

## Test the result when as.ratio is TRUE and Uniform distribution
test.syntheticNucMapFromDist_good_result_04 <- function() {
    obs <- syntheticNucMapFromDist(wp.num = 5, wp.del = 2,
                                   wp.var = 3, fuz.num = 2, fuz.var = 40,
                                   lin.len = 4, as.ratio = TRUE,
                                   rnd.seed = 26, distr = "Uniform")

    exp.wp.starts  <- c(1, 152, 454)
    exp.wp.nreads  <- c(3, 30, 80)
    exp.fuz.starts <- c(96, 326)
    exp.fuz.nreads <- c(91, 13)
    exp.nuc.len    <- 147

    message     <- paste0(" test.syntheticNucMapFromDist_good_result_04() ",
                          "- syntheticNucMapFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp.starts, exp.wp.starts, msg = message)
    checkEqualsNumeric(obs$wp.nreads, exp.wp.nreads, msg = message)
    checkEqualsNumeric(obs$fuz.starts, exp.fuz.starts, msg = message)
    checkEqualsNumeric(obs$fuz.nreads, exp.fuz.nreads, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$wp.reads), 113, msg = message)
    checkEqualsNumeric(length(obs$fuz.reads), 104, msg = message)
    checkEqualsNumeric(length(obs$syn.reads), 217, msg = message)
    checkEqualsNumeric(length(obs$syn.ratio), 664, msg = message)
    checkEqualsNumeric(length(obs$ctr.reads), 217, msg = message)
}

## Test the result when no nucleosome deleted
test.syntheticNucMapFromDist_good_result_05 <- function() {
    obs <- syntheticNucMapFromDist(wp.num = 5, wp.del = 0, max.cover = 15,
                                   wp.var = 23, fuz.num = 2, fuz.var = 43,
                                   lin.len = 12, as.ratio = FALSE,
                                   rnd.seed = 288, distr = "Student")

    exp.wp.starts  <- c(1, 160, 319, 478, 637)
    exp.wp.nreads  <- c(3, 11, 4, 11, 3)
    exp.fuz.starts <- c(660, 452)
    exp.fuz.nreads <- c(13, 8)
    exp.nuc.len    <- 147

    message     <- paste0(" test.syntheticNucMapFromDist_good_result_05() ",
                          "- syntheticNucMapFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp.starts, exp.wp.starts, msg = message)
    checkEqualsNumeric(obs$wp.nreads, exp.wp.nreads, msg = message)
    checkEqualsNumeric(obs$fuz.starts, exp.fuz.starts, msg = message)
    checkEqualsNumeric(obs$fuz.nreads, exp.fuz.nreads, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$wp.reads), 32, msg = message)
    checkEqualsNumeric(length(obs$fuz.reads), 21, msg = message)
    checkEqualsNumeric(length(obs$syn.reads), 53, msg = message)
    checkEqualsNumeric(length(obs$syn.ratio), 0, msg = message)
    checkEqualsNumeric(length(obs$ctr.reads), 0, msg = message)
}

