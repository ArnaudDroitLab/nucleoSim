#' @title Parameter validation for the \code{\link{syntheticNucMapFromDist}}
#' function
#'
#' @param wp.num the number of well-positioned
#' (non overlapping) nucleosomes.
#'
#' @param wp.del the number of well-positioned nucleosomes to remove to create
#' uncovered regions.
#'
#' @param wp.var the maximum variance in basepairs of the well-positioned
#' nucleosomes. This parameter introduces some variation in the position of
#' the reads describing a nucleosome.
#'
#' @param fuz.num the number of fuzzy nucleosomes. Those nucleosomes are
#' distributed accordingly to selected distribution all over the region.
#' Those nucleosomes can overlap other well-positioned or fuzzy nucleosomes.
#'
#' @param fuz.var the maximum variance of the fuzzy nucleosomes. This variance
#' can be different than the one used for the well-positioned
#' nucleosome reads.
#'
#' @param max.cover the maximum coverage for one nucleosome. The final coverage
#' can have a higher value than \code{max.cover} since reads from differents
#' nucleosomes can be overlapping.
#'
#' @param nuc.len the nucleosome length.
#'
#' @param lin.len the length of the DNA linker DNA.
#'
#' @param rnd.seed a single value, interpreted as an \code{integer}, or
#' \code{NULL}. If a \code{integer} is given, the value is used to set the seed
#' of the random number generator. By fixing the seed, the generated results
#' can be reproduced.
#'
#' @param as.ratio a \code{logical}, if \code{TRUE},  a synthetic naked DNA
#' control map is created and the ratio between it and the nucleosome coverage
#' are calculated. It can be used to simulate hybridization ratio data, like
#' the one in Tiling Arrays. Both control map and calculated ratio are
#' returned.
#'
#' @param show.plot a \code{logical}, if \code{TRUE}, a plot of the coverage
#' map, with the nucleosome calls and optionally the calculated ratio is
#' outputed.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @author Astrid Louise Deschenes
#' @keywords internal
syntheticNucMapFromDistValidation <- function(wp.num, wp.del, wp.var, fuz.num,
                                              fuz.var,
                                              max.cover,
                                              nuc.len,
                                              lin.len,
                                              rnd.seed,
                                              as.ratio,
                                              show.plot)
{

    if (!isInteger(wp.num) || wp.num < 0) {
        stop("wp.num must be a non-negative integer")
    }

    if (!isInteger(wp.del) || wp.del < 0) {
        stop("wp.del must be a non-negative integer")
    }

    if (!isInteger(wp.var) || wp.var < 0) {
        stop("wp.var must be a non-negative integer")
    }

    if (!isInteger(fuz.num) || fuz.num < 0) {
        stop("fuz.num must be a non-negative integer")
    }

    if (!isInteger(fuz.var) || fuz.var < 0) {
        stop("fuz.var must be a non-negative integer")
    }

    return(0)
}


#' @title Validate if a value is an integer
#'
#' @description Validate if the value passed to the function is an integer or
#' not. To be considered as an integer, the value must have a length
#' of 1. The type of value can be a \code{integer} or \code{numerical}.
#' However, a \code{numerical} must have the same value
#' once casted to a \code{integer}.  A \code{vector} of
#' integers will returned \code{FALSE}.
#'
#' @param value an object to validate.
#'
#' @return \code{TRUE} is the parameter is a integer; otherwise \code{FALSE}
#'
#' @author Astrid Louise Deschenes
#' @keywords internal
isInteger <- function(value) {
    return((is.integer(value) && length(value) == 1) || (is.numeric(value) &&
            as.integer(value) == value) && length(value) == 1)
}
