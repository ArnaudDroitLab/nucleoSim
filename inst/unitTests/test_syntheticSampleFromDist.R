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

## Test the result when a string is passed as offset parameter
test.syntheticSampleFromDist_string_offset <- function() {
    obs <- tryCatch(syntheticSampleFromDist(wp.num = 2, wp.del = 0,
                        wp.var = 30, fuz.num = 10, fuz.var = 40,
                        rnd.seed = 15, offset = "2",
                        distr = "Normal"), error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" syntheticSampleFromDist_string_offset() ",
                      "- A string as offset parameter did not generated ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}
