###################################################
# Created by Astrid Deschenes
# 2015-07-24
###################################################

###################################################
## Test the syntheticSampleFromDist function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "nucleoSim" )
}

### }}}

###########################################################
## Test the syntheticSampleFromDist() function parameters
###########################################################

################################
## wp.num
################################

## Test the result when a string is passed as wp.num parameter
test.syntheticSampleFromDist_string_wp_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = "2", wp.del = 0,
                                wp.var = 30, fuz.num = 10, fuz.var = 40,
                                rnd.seed = 15, offset = 2,
                                distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_wp_num() ",
                        "- A string as wp.num parameter did not generated ",
                        "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.num parameter
test.syntheticSampleFromDist_negative_wp_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = -1, wp.del = 0,
                                        wp.var = 30, fuz.num = 10, fuz.var = 40,
                                        rnd.seed = 15, offset = 2,
                                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_wp_num() ",
                      "- A negative integer as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.num parameter
test.syntheticSampleFromDist_vector_wp_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = c(1, 2), wp.del = 0,
                                            wp.var = 30, fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 3,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticSampleFromDist_vector_wp_num() ",
                      "- A vector as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## wp.del
################################

## Test the result when a string is passed as wp.del parameter
test.syntheticSampleFromDist_string_wp_del <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = "0",
                                wp.var = 30, fuz.num = 10, fuz.var = 40,
                                rnd.seed = 15, offset = 2,
                                distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_wp_del() ",
                      "- A string as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.del parameter
test.syntheticSampleFromDist_negative_wp_del <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 12, wp.del = -1,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_wp_del() ",
                      "- A negative integer as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.del parameter
test.syntheticSampleFromDist_vector_wp_del <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = c(1, 4),
                                            wp.var = 30, fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_wp_del() ",
                      "- A vector as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## wp.var
################################

## Test the result when a string is passed as wp.var parameter
test.syntheticSampleFromDist_string_wp_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                                    wp.var = "30", fuz.num = 10, fuz.var = 40,
                                    rnd.seed = 15, offset = 2,
                                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_wp_var() ",
                      "- A string as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as wp.var parameter
test.syntheticSampleFromDist_negative_wp_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = -1, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_wp_var() ",
                      "- A negative integer as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as wp.var parameter
test.syntheticSampleFromDist_vector_wp_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = c(1, 3), fuz.num = 3, fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_wp_var() ",
                      "- A vector as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## fuz.num
################################

## Test the result when a string is passed as fuz.num parameter
test.syntheticSampleFromDist_string_fuz_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                                    wp.var = 30, fuz.num = "10", fuz.var = 40,
                                    rnd.seed = 15, offset = 2,
                                    distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_fuz_num() ",
                      "- A string as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.num parameter
test.syntheticSampleFromDist_negative_fuz_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = -1, fuz.var = 40,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_fuz_num() ",
                      "- A negative integer as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.num parameter
test.syntheticSampleFromDist_vector_fuz_num <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = c(1, 3), fuz.var = 40,
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_fuz_num() ",
                      "- A vector as fuz.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## fuz.var
################################

## Test the result when a string is passed as fuz.var parameter
test.syntheticSampleFromDist_string_fuz_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = "40",
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_fuz_var() ",
                      "- A string as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as fuz.var parameter
test.syntheticSampleFromDist_negative_fuz_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = 12, fuz.var = -1,
                                            rnd.seed = 15, offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_fuz_var() ",
                      "- A negative integer as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as fuz.var parameter
test.syntheticSampleFromDist_vector_fuz_var <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = c(1, 3),
                                            rnd.seed = 15, offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_fuz_var() ",
                      "- A vector as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## rnd.seed
################################

## Test the result when a string is passed as rnd.seed parameter
test.syntheticSampleFromDist_string_rnd_seed <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = "15", offset = 2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticSampleFromDist_string_rnd_seed() ",
                      "- A string as rnd.seed parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as rnd.seed parameter
test.syntheticSampleFromDist_vector_rnd_seed <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = 2,
                                            rnd.seed = c(1, 2), offset = 12,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "rnd.seed must be NULL or an integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_rnd_seed() ",
                      "- A vector as rnd.seed parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## offset
################################

## Test the result when a string is passed as offset parameter
test.syntheticSampleFromDist_string_offset <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                                            rnd.seed = 15, offset = "2",
                                            distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_string_offset() ",
                      "- A string as offset parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as offset parameter
test.syntheticSampleFromDist_negative_offset <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 12, wp.del = 1,
                                            wp.var = 1, fuz.num = 12, fuz.var = 1,
                                            rnd.seed = 15, offset = -2,
                                            distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_negative_offset() ",
                      "- A negative integer as offset  parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a vector is passed as offset parameter
test.syntheticSampleFromDist_vector_offset <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 11, wp.del = 1,
                                            wp.var = 4, fuz.num = 1, fuz.var = 1,
                                            rnd.seed = 15, offset = c(1,2),
                                            distr = "Student"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticSampleFromDist_vector_offset() ",
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

## Test the result when as.ratio is FALSE
test.syntheticSampleFromDist_good_result_01 <- function() {
    obs <- syntheticSampleFromDist(wp.num = 4, wp.del = 2,
                                    wp.var = 3, fuz.num = 1, fuz.var = 40,
                                    lin.len = 4, rnd.seed = 125, offset = 100,
                                    distr = "Normal")

    exp.wp <- data.frame(nucleopos=c(175, 326, 477, 628), nreads=c(83, 13, 31, 0))
    exp.fuz <- data.frame(nucleopos=c(540), nreads=c(48))

    message     <- paste0(" test.syntheticSampleFromDist_good_result_01() ",
                          "- syntheticSampleFromDist did not generated ",
                          "expected values")

    checkEqualsNumeric(obs$wp, exp.wp, msg = message)
    checkEqualsNumeric(obs$fuz, exp.fuz, msg = message)
    checkEqualsNumeric(length(obs$dataIP), 4, msg = message)
    checkEqualsNumeric(nrow(obs$dataIP), 350, msg = message)
    checkEqualsNumeric(length(obs$paired), 3, msg = message)
    checkEqualsNumeric(nrow(obs$paired), 175, msg = message)
}
