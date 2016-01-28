###################################################
# Created by Astrid Deschenes
# 2015-07-24
###################################################

###################################################
## Test the syntheticNucReadsFromDist function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "nucleoSim" )
}

### }}}

###########################################################
## Test the syntheticNucReadsFromDist() function parameters
###########################################################

################################
## wp.num
################################

## Test the result when a string is passed as wp.num parameter
test.syntheticNucReadsFromDist_string_wp_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = "2", wp.del = 0,
                                wp.var = 30, fuz.num = 10, fuz.var = 40,
                                rnd.seed = 15, offset = 2,
                                distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_num() ",
                        "- A string as wp.num parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.num parameter
test.syntheticNucReadsFromDist_negative_wp_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = -1, wp.del = 0,
                                        wp.var = 30, fuz.num = 10, fuz.var = 40,
                                        rnd.seed = 15, offset = 2,
                                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_num() ",
                      "- A negative integer as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.num parameter
test.syntheticNucReadsFromDist_vector_wp_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = c(1, 2), wp.del = 0,
                                            wp.var = 30, fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 3,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucReadsFromDist_vector_wp_num() ",
                      "- A vector as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## wp.del
################################

## Test the result when a string is passed as wp.del parameter
test.syntheticNucReadsFromDist_string_wp_del <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = "0",
                                wp.var = 30, fuz.num = 10, fuz.var = 40,
                                rnd.seed = 15, offset = 2,
                                distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_del() ",
                      "- A string as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.del parameter
test.syntheticNucReadsFromDist_negative_wp_del <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = -1,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_del() ",
                      "- A negative integer as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.del parameter
test.syntheticNucReadsFromDist_vector_wp_del <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = c(1, 4),
                                            wp.var = 30, fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_wp_del() ",
                      "- A vector as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## wp.var
################################

## Test the result when a string is passed as wp.var parameter
test.syntheticNucReadsFromDist_string_wp_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                    wp.var = "30", fuz.num = 10, fuz.var = 40,
                                    rnd.seed = 15, offset = 2,
                                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_wp_var() ",
                      "- A string as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.var parameter
test.syntheticNucReadsFromDist_negative_wp_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = -1, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_wp_var() ",
                      "- A negative integer as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.var parameter
test.syntheticNucReadsFromDist_vector_wp_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = c(1, 3), fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_wp_var() ",
                      "- A vector as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## fuz.num
################################

## Test the result when a string is passed as fuz.num parameter
test.syntheticNucReadsFromDist_string_fuz_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                    wp.var = 30, fuz.num = "10", fuz.var = 40,
                                    rnd.seed = 15, offset = 2,
                                    distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_fuz_num() ",
                      "- A string as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.num parameter
test.syntheticNucReadsFromDist_negative_fuz_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = -1, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_fuz_num() ",
                      "- A negative integer as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.num parameter
test.syntheticNucReadsFromDist_vector_fuz_num <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = c(1, 3), fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_fuz_num() ",
                      "- A vector as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## fuz.var
################################

## Test the result when a string is passed as fuz.var parameter
test.syntheticNucReadsFromDist_string_fuz_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = "40",
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_fuz_var() ",
                      "- A string as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.var parameter
test.syntheticNucReadsFromDist_negative_fuz_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = 12, fuz.var = -1,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_fuz_var() ",
                      "- A negative integer as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.var parameter
test.syntheticNucReadsFromDist_vector_fuz_var <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = c(1, 3),
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_fuz_var() ",
                      "- A vector as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## rnd.seed
################################

## Test the result when a string is passed as rnd.seed parameter
test.syntheticNucReadsFromDist_string_rnd_seed <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = "15", offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_rnd_seed() ",
                      "- A string as rnd.seed parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as rnd.seed parameter
test.syntheticNucReadsFromDist_vector_rnd_seed <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = 2,
                                            rnd.seed = c(1, 2), offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_rnd_seed() ",
                      "- A vector as rnd.seed parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## offset
################################

## Test the result when a string is passed as offset parameter
test.syntheticNucReadsFromDist_string_offset <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = "2",
                                            distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_string_offset() ",
                      "- A string as offset parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as offset parameter
test.syntheticNucReadsFromDist_negative_offset <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = 12, fuz.var = 1,
                                            rnd.seed = 15, offset = -2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_negative_offset() ",
                      "- A negative integer as offset  parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as offset parameter
test.syntheticNucReadsFromDist_vector_offset <- function() {
    obs <- tryCatch(syntheticNucReadsFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = 1,
                                            rnd.seed = 15, offset = c(1,2),
                                            distr = "Student"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromDist_vector_offset() ",
                      "- A vector as offset parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## distr
################################



################################
## Good results
################################

## Test the result for a specific case
test.syntheticNucReadsFromDist_good_result_01 <- function() {
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 2,
                                    wp.var = 3, fuz.num = 1, fuz.var = 40,
                                    lin.len = 4, rnd.seed = 125, offset = 100,
                                    distr = "Normal")

    exp.wp <- data.frame(nucleopos=c(175, 326), nreads=c(83, 13))
    exp.fuz      <- data.frame(nucleopos=c(236), nreads=c(61))
    exp.nuc.len <- 147

    message     <- paste0(" test.syntheticNucReadsFromDist_good_result_01() ",
                          "- syntheticNucReadsFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp, exp.wp, msg = message)
    checkEqualsNumeric(obs$fuz, exp.fuz, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$dataIP), 4, msg = message)
    checkEqualsNumeric(nrow(obs$dataIP), 314, msg = message)
    checkEqualsNumeric(length(obs$paired), 3, msg = message)
    checkEqualsNumeric(nrow(obs$paired), 157, msg = message)
}

## Test the result for a specific case
test.syntheticNucReadsFromDist_good_result_02 <- function() {
    obs <- syntheticNucReadsFromDist(wp.num = 4, wp.del = 3, nuc.len = 144,
                                     wp.var = 6, fuz.num = 2, fuz.var = 60,
                                     lin.len = 4, rnd.seed = 129, offset = 1000,
                                     distr = "Student")

    exp.wp <- data.frame(nucleopos=c(1369), nreads=c(17))
    exp.fuz      <- data.frame(nucleopos=c(1357, 1418), nreads=c(98, 54))
    exp.nuc.len <- 144

    message     <- paste0(" test.syntheticNucReadsFromDist_good_result_02() ",
                          "- syntheticNucReadsFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp, exp.wp, msg = message)
    checkEqualsNumeric(obs$fuz, exp.fuz, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$dataIP), 4, msg = message)
    checkEqualsNumeric(nrow(obs$dataIP), 338, msg = message)
    checkEqualsNumeric(length(obs$paired), 3, msg = message)
    checkEqualsNumeric(nrow(obs$paired), 169, msg = message)
}
