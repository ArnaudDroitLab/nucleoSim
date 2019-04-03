###################################################
# Created by Astrid Deschenes
# 2016-02-04
###################################################

###################################################
## Test the syntheticNucReadsFromMap function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "nucleoSim" )
}

### }}}

NUCLEO_MAP <- syntheticNucMapFromDist(wp.num=3, wp.del=1, wp.var=12,
                                         fuz.num=1, fuz.var=44, max.cover=65,
                                         nuc.len=147, len.var=2, lin.len=40,
                                         rnd.seed=155, distr="Uniform")


###########################################################
## Test the syntheticNucReadsFromMap() function parameters
###########################################################

################################
## syntheticNucMap
################################

## Test the result when a string is passed as syntheticNucMap parameter
test.syntheticNucReadsFromMap_string_syntheticNucMap <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap="canada",
                            read.len=40, offset=12), error = conditionMessage)
    exp <- "syntheticNucMap must be an object of class \'syntheticNucMap\'"
    message <- paste0(" test.syntheticNucReadsFromMap_string_syntheticNucMap() ",
                      "- A string as syntheticNucMap parameter did not ",
                      "generate expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when an integer is passed as syntheticNucMap parameter
test.syntheticNucReadsFromMap_integer_syntheticNucMap <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=22,
                                             read.len=40, offset=12), error = conditionMessage)
    exp <- "syntheticNucMap must be an object of class \'syntheticNucMap\'"
    message <- paste0(" test.syntheticNucReadsFromMap_integer_syntheticNucMap() ",
                      "- An integer as syntheticNucMap parameter did not ",
                      "generate expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## read.len
################################

## Test the result when a string is passed as read.len parameter
test.syntheticNucReadsFromMap_string_read_len <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=NUCLEO_MAP,
                                                read.len="igloo", offset=12),
                                                error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromMap_string_read_len() ",
                      "- A string as read.len parameter did not generate ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative integer is passed as syntheticNucMap parameter
test.syntheticNucReadsFromMap_negative_integer_syntheticNucMap <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=NUCLEO_MAP,
                                                read.len=-1, offset=12),
                                                error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromMap_negative_integer_syntheticNucMap() ",
                      "- An negative integer as syntheticNucMap parameter ",
                      "did not generate expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zero is passed as syntheticNucMap parameter
test.syntheticNucReadsFromMap_zero_syntheticNucMap <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=NUCLEO_MAP,
                                                read.len=0, offset=12),
                                                error = conditionMessage)
    exp <- "read.len must be a positive integer"
    message <- paste0(" test.syntheticNucReadsFromMap_zero_syntheticNucMap() ",
                      "- A zero as syntheticNucMap parameter ",
                      "did not generate expected error.")
    checkEquals(obs, exp, msg = message)
}

################################
## offset
################################

## Test the result when a string is passed as read.len parameter
test.syntheticNucReadsFromMap_string_offset <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=NUCLEO_MAP,
                                                read.len=30, offset="testing"),
                                                error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromMap_string_offset() ",
                      "- A string as offset parameter did not generate ",
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a negative integer is passed as offset parameter
test.syntheticNucReadsFromMap_negative_integer_offset <- function() {
    obs <- tryCatch(syntheticNucReadsFromMap(syntheticNucMap=NUCLEO_MAP,
                                                read.len=20, offset=-1),
                                                error = conditionMessage)
    exp <- "offset must be a non-negative integer"
    message <- paste0(" test.syntheticNucReadsFromMap_negative_integer_offset() ",
                      "- An negative integer as offset parameter ",
                      "did not generate expected error.")
    checkEquals(obs, exp, msg = message)
}


################################
## Good results
################################

## Test the result for a specific case
test.syntheticNucReadsFromMap_good_result_01 <- function() {
    obs <- syntheticNucReadsFromMap(syntheticNucMap = NUCLEO_MAP,
                                                        read.len = 40,
                                                        offset = 1000)

    exp.wp <- data.frame(nucleopos=c(1262, 1449), nreads=c(29, 50))
    exp.fuz <- data.frame(nucleopos=c(1481), nreads=c(25))
    exp.nuc.len <- 147
    exp.dataIP.colnames <- c("chr", "start", "end", "strand", "ID")
    exp.dataIP.chr <- rep("chr_SYNTHETIC", 208)
    exp.paired.colnames <- c("chr", "start", "end", "ID")
    exp.paired.chr <- rep("chr_SYNTHETIC", 104)

    message     <- paste0(" test.syntheticNucReadsFromMap_good_result_01() ",
                          "- syntheticNucReadsFromMap did not generate ",
                          "expected values")

    checkEqualsNumeric(obs$wp, exp.wp, msg = message)
    checkEqualsNumeric(obs$fuz, exp.fuz, msg = message)
    checkEqualsNumeric(obs$nuc.len, exp.nuc.len, msg = message)
    checkEqualsNumeric(length(obs$dataIP), 5, msg = message)
    checkEqualsNumeric(colnames(obs$dataIP), exp.dataIP.colnames, msg = message)
    checkEqualsNumeric(nrow(obs$dataIP), 208, msg = message)
    checkEquals(as.character(obs$dataIP$chr), exp.dataIP.chr, msg = message)
    checkEqualsNumeric(length(obs$paired), 4, msg = message)
    checkEqualsNumeric(colnames(obs$paired), exp.paired.colnames, msg = message)
    checkEqualsNumeric(nrow(obs$paired), 104, msg = message)
    checkEquals(as.character(obs$paired$chr), exp.paired.chr, msg = message)

}

