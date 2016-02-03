#' @title Parameter validation for the \code{\link{syntheticNucMapFromDist}}
#' function
#'
#' @description Validate that all values passed to the function are formated
#' for the \code{\link{syntheticNucMapFromDist}} function.
#'
#' @param wp.num a non-negative \code{integer}, the number of well-positioned
#' (non overlapping) nucleosomes.
#'
#' @param wp.del a non-negative \code{integer}, the number of well-positioned
#' nucleosomes to remove to create uncovered regions.
#'
#' @param wp.var a non-negative \code{integer}, the variance in base pairs
#' of the well-positioned
#' nucleosomes. This parameter introduces some variation in the position of
#' the reads describing a nucleosome.
#'
#' @param fuz.num a non-negative \code{integer}, the number of fuzzy
#' nucleosomes. Those nucleosomes are distributed accordingly to selected
#' distribution all over
#' the region. Those nucleosomes can overlap other well-positioned or
#' fuzzy nucleosomes.
#'
#' @param fuz.var a non-negative \code{integer}, the variance of the fuzzy
#' nucleosomes. This variance
#' can be different than the one used for the well-positioned
#' nucleosome reads.
#'
#' @param max.cover a positive \code{integer}, the maximum coverage for one
#' nucleosome. The final coverage
#' can have a higher value than \code{max.cover} since reads from differents
#' nucleosomes can be overlapping.
#'
#' @param nuc.len a non-negative \code{integer}, the nucleosome length.
#'
#' @param len.var a non-negative \code{integer}, the variance of the
#' distance between a forward read and its paired reverse read.
#'
#' @param lin.len a non-negative \code{integer}, the length of the DNA
#' linker DNA.
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
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @author Astrid Deschenes
#' @keywords internal
syntheticNucMapFromDistValidation <- function(wp.num, wp.del, wp.var, fuz.num,
                                                fuz.var,
                                                max.cover,
                                                nuc.len,
                                                len.var,
                                                lin.len,
                                                rnd.seed,
                                                as.ratio)
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

    if (!isInteger(max.cover) || max.cover < 1) {
        stop("max.cover must be a positive integer")
    }

    if (!isInteger(nuc.len) || nuc.len < 1) {
        stop("nuc.len must be a positive integer")
    }

    if (!isInteger(len.var) || len.var < 0) {
        stop("len.var must be a non-negative integer")
    }

    if (!isInteger(lin.len) || lin.len < 0) {
        stop("lin.len must be a non-negative integer")
    }

    if (!is.null(rnd.seed) && !isInteger(rnd.seed)) {
        stop("rnd.seed must be NULL or an integer")
    }

    if (!is.logical(as.ratio)) {
        stop("as.ratio must be a logical value")
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
#' @author Astrid Deschenes
#' @keywords internal
isInteger <- function(value) {
    return((is.integer(value) && length(value) == 1) || (is.numeric(value) &&
            as.integer(value) == value) && length(value) == 1)
}


#' @title Subsection of parameter validation for identical parameters between
#' \code{syntheticNucReadsFromMap} and \code{syntheticNucReadsFromDist}
#' functions.
#'
#' @description Validate that identical values passed to both
#' \code{syntheticNucReadsFromMap} and \code{syntheticNucReadsFromDist}
#' functions are correctly formatted.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @author Astrid Deschenes
#' @keywords internal
syntheticNucReadsValidation <- function(read.len, offset)
{
    ## Validate that offset is a positive integer
    if (!isInteger(read.len) || read.len <= 0) {
        stop("read.len must be a positive integer")
    }

    ## Validate that offset is a non-negative integer
    if (!isInteger(offset) || offset < 0) {
        stop("offset must be a non-negative integer")
    }

    return(0)
}


#' @title Create a synthetic nucleosome reads from a synthetic nucleosome map
#'
#' @description Generate a synthetic nucleosome map using a synthetic
#' nucleosome map.
#'
#' @param map a \code{list} of \code{class} "syntheticNucMap"
#'
#' @param read.len the length of each of the paired-end reads. Default = 40.
#'
#' @param offset the number of bases used to offset all nucleosomes and reads.
#' This is done to ensure that all nucleosome positions and read alignment
#' are of positive values.
#'
#' @param call the function call.
#'
#' @return an \code{list} of \code{class} "syntheticNucReads" containing the
#' following elements:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{dataIP} a \code{data.frame} with the chromosome name, the
#' starting and ending positions and the direction of all forward
#' and reverse reads for all well-positioned and fuzzy nucleosomes.
#' \item \code{wp} a \code{data.frame} with the positions of all the
#' well-positioned nucleosomes, as well as the number of paired-reads
#' associated to each one.
#' \item \code{fuz} a \code{data.frame} with the positions of all the fuzzy
#' nucleosomes, as well as the number of paired-reads associated to each one.
#' \item \code{paired} a \code{data.frame} with the starting and ending
#' positions of the reads used to generate the paired-end reads.
#' }
#'
#' @author Astrid Deschenes
#' @keywords internal
createNucReadsFromNucMap <- function(map, read.len, offset, call) {


    ## Extract reads to create a dataframe
    syn.reads <- as.data.frame(map$syn.reads)

    ## Get number of reads present in dataframe
    nreads <- nrow(syn.reads)

    ## Add an offset to the read positions
    syn.reads$start <- syn.reads$start + offset
    syn.reads$end   <- syn.reads$end   + offset

    # Store read paired information
    paired <- data.frame(rep("chr1", nreads),
                         as.integer(syn.reads$start),
                         as.integer(syn.reads$end))
    colnames(paired) <- c("chr", "start", "end")

    # Order paired reads by starting position
    paired <- paired[order(paired$start), ]

    ## Create dataset with forward and reverse reads for all nucleosomes
    ## Forward and reverse reads are generated from read paired information
    ## using the read.len parameter
    dataIP <- data.frame(rep("chr_SYNTHETIC", 2 * nreads),
                         c(as.integer(syn.reads$start),
                           as.integer(syn.reads$end - read.len)),
                         c(as.integer(syn.reads$start + read.len),
                           as.integer(syn.reads$end)),
                         c(rep("+", nreads), rep("-", nreads)))
    colnames(dataIP) = c("chr", "start", "end", "strand")

    # Order reads by starting position
    dataIP <- dataIP[order(dataIP$start), ]
    rownames(dataIP)<-1:nrow(dataIP)

    # Create data.frame with well-positioned nucleosomes shifted by offset value
    space <- round(map$nuc.len/2)
    wp <- data.frame(map$wp.starts + offset + space, map$wp.nreads)
    colnames(wp) <- c("nucleopos", "nreads")

    # Create data.frame with fuzzy nucleosomes shifted by offset value
    fuz <- data.frame(map$fuz.starts + offset + space,
                      map$fuz.nreads)
    colnames(fuz) <- c("nucleopos", "nreads")

    # Create returned result list
    result <- list(call = call, dataIP = dataIP, wp = wp, fuz = fuz,
                   paired = paired, nuc.len = map$nuc.len)

    # Return a list marked as an syntheticNucReads class
    class(result) <- "syntheticNucReads"
    return(result)
}
