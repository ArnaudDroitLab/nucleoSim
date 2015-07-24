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

