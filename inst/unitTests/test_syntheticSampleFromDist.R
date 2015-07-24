###################################################
# Created by Astrid Louise Deschenes
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

## Test the result when a string is passed as wp.num parameter
test.syntheticNucMapFromDist_string_wp_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = "allo", wp.del = 0,
                            wp.var = 30, fuz.num = 10, fuz.var = 40,
                            as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "wp.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_fuz_num() ",
                      "- A vector as wp.num parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a string is passed as fuz.num parameter
test.syntheticNucMapFromDist_string_fuz_num <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 0,
                        wp.var = 30, fuz.num = "allo", fuz.var = 40,
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                            as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                distr = "Normal"), error = conditionMessage)
    exp <- "fuz.num must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_fuz_num() ",
                "- A vector as fuz.num parameter did not generated ",
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a string is passed as wp.del parameter
test.syntheticNucMapFromDist_string_wp_del <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = "allo",
                    wp.var = 30, fuz.num = 4, fuz.var = 40,
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.del must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_wp_del() ",
                      "- A vector as wp.del parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

