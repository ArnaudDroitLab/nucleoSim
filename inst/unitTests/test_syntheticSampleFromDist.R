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

## Test the result when a string is passed as wp.var parameter
test.syntheticNucMapFromDist_string_wp_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                            wp.var = "allo", fuz.num = 4, fuz.var = 40,
                            as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "wp.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_wp_var() ",
                      "- A vector as wp.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a string is passed as fuz.var parameter
test.syntheticNucMapFromDist_string_fuz_var <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                    wp.var = 2, fuz.num = 4, fuz.var = "allo",
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                            as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                            distr = "Normal"), error = conditionMessage)
    exp <- "fuz.var must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_fuz_var() ",
                      "- A vector as fuz.var parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a string is passed as max.cover parameter
test.syntheticNucMapFromDist_string_max_cover <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 3, wp.del = 2,
                wp.var = 2, fuz.num = 4, fuz.var = 2, max.cover = "allo",
                as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
            as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                        distr = "Normal"), error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_string_lin_len() ",
                      "- A string as lin.len parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative is passed as lin.len parameter
test.syntheticNucMapFromDist_negative_lin_len <- function() {
    obs <- tryCatch(syntheticNucMapFromDist(wp.num = 10, wp.del = 1,
                        wp.var = 1, fuz.num = 4, fuz.var = 40, lin.len = -1,
                        as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
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
                    as.ratio = TRUE, show.plot = FALSE, rnd.seed = 15,
                    distr = "Normal"), error = conditionMessage)
    exp <- "lin.len must be a non-negative integer"
    message <- paste0(" syntheticNucMapFromDist_vector_lin_len() ",
                      "- A vector as lin.len parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

