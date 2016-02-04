#' @title Generate a synthetic nucleosome map
#'
#' @description Generate a synthetic nucleosome map using the distribution
#' selected by the user. The user has choice between three different
#' distributions: Normal, Student and Uniform.
#'
#' This function is a modified version of the syntheticNucMap() function
#' from Bioconductor nucleR package (Flores and Orozco, 2011).
#'
#' @param wp.num a non-negative \code{integer}, the number of well-positioned
#' (non overlapping) nucleosomes.
#'
#' @param wp.del a non-negative \code{integer}, the number of well-positioned
#' nucleosomes to remove to create uncovered regions.
#'
#' @param wp.var a non-negative \code{integer}, the maximum variance in
#' basepairs of the well-positioned
#' nucleosomes. This parameter introduces some variation in the position of
#' the reads describing a nucleosome.
#'
#' @param fuz.num a non-negative \code{integer}, the number of
#' fuzzy nucleosomes. Those nucleosomes are
#' distributed accordingly to selected distribution all over the region.
#' Those nucleosomes can overlap other well-positioned or fuzzy nucleosomes.
#'
#' @param fuz.var a non-negative \code{integer}, the maximum variance of the
#' fuzzy nucleosomes. This variance
#' can be different than the one used for the well-positioned
#' nucleosome reads.
#'
#' @param max.cover a positive \code{integer}, the maximum coverage for
#' one nucleosome. The final coverage
#' can have a higher value than \code{max.cover} since reads from differents
#' nucleosomes can be overlapping. Default = 100.
#'
#' @param nuc.len a non-negative \code{numeric}, the nucleosome length.
#' Default = 147.
#'
#' @param len.var a non-negative \code{integer}, the variance of the
#' distance between a forward read and its paired reverse read. Default = 10.
#'
#' @param lin.len a non-negative \code{integer}, the length of
#' the DNA linker. Default = 20.
#'
#' @param rnd.seed a single value, interpreted as an \code{integer}, or
#' \code{NULL}. If a \code{integer} is given, the value is used to set the seed
#' of the random number generator. By fixing the seed, the generated results
#' can be reproduced. Default = \code{NULL}.
#'
#' @param as.ratio a \code{logical}, if \code{TRUE}, a synthetic naked DNA
#' control map is created and the ratio between it and the nucleosome coverage
#' are calculated. It can be used to simulate hybridization ratio data, like
#' the one in Tiling Arrays. Both control map and calculated ratio are
#' returned. Default = \code{FALSE}.
#'
#' @param distr the name of the distribution used to generate the nucleosome
#' map. The choices are : \code{"Uniform"}, \code{"Normal"} and
#' \code{"Student"}. Default = \code{"Uniform"}.
#'
#' @return an \code{list} of \code{class} "syntheticNucMap" containing the
#' following elements:
#' \itemize{
#'     \item \code{call} the matched call.
#'     \item \code{wp.starts} a \code{vector} of \code{integer}, the start
#' positions of all well-positioned nucleosome regions. The central
#' position of the nucleosome is calculated as wp.starts + round(nuc.len/2).
#'     \item \code{wp.nreads} a \code{vector} of \code{integer}, the number of
#' reads associated to each well-positioned nucelosome.
#'     \item \code{wp.reads} a \code{IRanges} containing the well positioned
#' nucleosome reads.
#'     \item \code{fuz.starts} a \code{vector} of \code{integer}, the
#' start position of all the fuzzy nucleosomes.
#'     \item \code{fuz.nreads} a \code{vector} of \code{integer}, the number
#' of reads associated to each fuzzy nucleosome.
#'     \item \code{fuz.reads} a \code{IRanges} containing the fuzzy nucleosome
#' reads.
#'     \item \code{syn.reads} a \code{IRanges} containing all the synthetic
#' nucleosome reads (from both fuzzy and well-positioned nucleosomes).
#'     \item \code{nuc.len} a \code{numeric} the nucleosome length.
#' }
#' The following elements will be only returned if \code{as.ratio=TRUE}:
#' \itemize{
#'     \item \code{ctr.reads} a \code{IRanges} containing the naked DNA
#' (control) reads.
#'     \item \code{syn.ratio} a \code{Rle} containing the calculated ratio
#' between the nucleosome coverage and the control coverage.
#' }
#'
#' @examples
#'
#' ## Generate a synthetic map with 20 well-positioned nucleosomes and 10 fuzzy
#' ## nucleosomes using a Normal distribution with a variance of 30 for the
#' ## well-positioned nucleosomes, a variance of 40 for the fuzzy nucleosomes
#' ## and a seed of 15
#' syntheticNucMapFromDist(wp.num = 20, wp.del = 0, wp.var = 30,
#'     fuz.num = 10, fuz.var = 40, rnd.seed = 15,
#'     distr = "Normal")
#'
#' ## Same output but with ratio
#' syntheticNucMapFromDist(wp.num = 20, wp.del = 0, wp.var = 30,
#'     fuz.num = 10, fuz.var = 40,
#'     rnd.seed = 15, as.ratio = TRUE, distr = "Normal")
#'
#'
#' @author Rawane Samb, Astrid Deschenes
#' @importFrom stats runif rnorm rt
#' @importFrom IRanges IRanges coverage start
#' @importFrom graphics plot lines abline points legend
#' @importFrom S4Vectors Rle start
#' @export
syntheticNucMapFromDist <- function(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                                max.cover = 100,
                                nuc.len = 147,
                                len.var = 10,
                                lin.len = 20,
                                rnd.seed = NULL,
                                as.ratio = FALSE,
                                distr= c("Uniform", "Normal", "Student"))
{
    # Get call information
    cl <- match.call()

    # Validate parameters
    syntheticNucMapFromDistValidation(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                            max.cover, nuc.len, len.var, lin.len, rnd.seed,
                            as.ratio)

    # Validate distribution value
    distr <- match.arg(distr, several.ok = FALSE)

    # Set seed if given
    if (!is.null(rnd.seed)) {
        set.seed(rnd.seed)
    }

    ###########################################################################
    # WELL POSITIONED NUCLEOSOMES
    ###########################################################################

    # Starting point of putative nucleosomes
    wp.starts <- (nuc.len + lin.len) * seq(0, wp.num - 1) + 1

    # How many times a read is repeated
    wp.nreads <- round(runif(wp.num, min = 1, max = max.cover))

    # Delete some nucleosomes (set repetition time to 0)
    wp.deleted <- sample(x = 1:wp.num, size = wp.del, replace = FALSE)
    wp.nreads[wp.deleted] <- 0

    # Set each nucleosome as a repeated single start position
    wp.repstar <- rep(wp.starts, wp.nreads)

    # Add some variance to the starting points
    if (distr == "Normal")
    {
        wp.varstar <- wp.repstar + round(rnorm(length(wp.repstar), 0,
                                                wp.var^0.5))
    }
    else if (distr == "Student")
    {
        wp.varstar <- wp.repstar + round(rt(length(wp.repstar), 4))
    }
    else
    {
        wp.varstar <- wp.repstar + round(runif(length(wp.repstar),
                                                min=-wp.var, max=wp.var))
    }

    # Add some variance to the length of the sequences
    # The variance associated to the length of the sequences is always
    # following a normal distribution
    wp.varlen <- nuc.len + round(rnorm(length(wp.repstar), 0, len.var^0.5))

    # Putative reads
    wp.reads <- IRanges(start=wp.varstar, end=wp.varstar + wp.varlen)

    ###########################################################################
    # OVERLAPPED (FUZZY) NUCLEOSOMES
    ###########################################################################

    # Starting point of fuzzy nucleosomes (random)
    fuz.starts <- round(runif(fuz.num, min = 1,
                                    max = (nuc.len + lin.len) * wp.num))

    # How many times a read is repeated
    fuz.nreads <- round(runif(fuz.num, min = 1, max = max.cover))

    # Set each nucleosome as a repeated single start position
    fuz.repstar <- rep(fuz.starts, fuz.nreads)

    # Add some variance to the starting points
    if (distr == "Normal")
    {
        fuz.varstar <- fuz.repstar + round(rnorm(length(fuz.repstar), 0,
                                                fuz.var^0.5))
    }
    else if (distr == "Student")
    {
        fuz.varstar <- fuz.repstar + round(rt(length(fuz.repstar), 4))
    }
    else
    {
        fuz.varstar <- fuz.repstar + round(runif(length(fuz.repstar),
                                                min = -fuz.var,
                                                max = fuz.var))
    }

    # Add some variance to the length of the sequences
    # The variance associated to the length of the sequences is always
    # following a normal distribution
    fuz.varlen  <- nuc.len + round(rnorm(length(fuz.repstar), 0, len.var^0.5))

    # Overlapped reads
    fuz.reads <- IRanges(start = fuz.varstar, end = fuz.varstar + fuz.varlen)

    # All synthetic reads (from fuzzy and well-positioned nucleosomes)
    syn.reads <- c(wp.reads, fuz.reads)

    # RATIO AS HYBRIDIZATION (Tiling Array)
    if (as.ratio) {
        # Just put the same amount of reads as before randomly
        ctr.starts <- round(runif(length(syn.reads),
                                    min = 1,
                                    max = max(start(syn.reads))))

        # A random read length, between 50 and 250
        ctr.widths <- round(runif(length(syn.reads), min = 50, max = 250))

        # Create control reads
        ctr.reads <- IRanges(start = ctr.starts, width = ctr.widths)

        # Calculate ratio
        syn.ratio <- suppressWarnings(log2(as.vector(coverage(syn.reads))) -
                                        log2(as.vector(coverage(ctr.reads))))

        # Some lost bases... such as in reality
        syn.ratio[abs(syn.ratio) == Inf] <- NA
        syn.ratio <- Rle(syn.ratio)
    }

    # Preparing returned list
    result <- list(call = cl)

    # The deleted nucleosomes are not returned
    if (wp.del > 0) {
        result[["wp.starts"]] <- wp.starts[-c(wp.deleted)]
        result[["wp.nreads"]] <- wp.nreads[-c(wp.deleted)]
    } else {
        result[["wp.starts"]] <- wp.starts
        result[["wp.nreads"]] <- wp.nreads
    }
    result[["wp.reads"]]  <- wp.reads

    result[["fuz.starts"]] <- fuz.starts
    result[["fuz.nreads"]] <- fuz.nreads
    result[["fuz.reads"]]  <- fuz.reads

    result[["syn.reads"]]  <- syn.reads

    if(as.ratio) {
        result[["ctr.reads"]] <- ctr.reads
        result[["syn.ratio"]] <- syn.ratio
    }

    result[["nuc.len"]] <- nuc.len

    # Return a list marked as an syntheticNucMap class
    class(result) <- "syntheticNucMap"
    return(result)
}


#' @title Generate a synthetic nucleosome map containing forward and
#' reverse reads
#'
#' @description Generate a synthetic nucleosome map using the distribution
#' selected by the user. The user has choice between three different
#' distributions: Normal, Student and Uniform. The final map is composed of
#' paired-end reads.
#'
#' This function is using a modified version of the syntheticNucMap() function
#' from Bioconductor nucleR package (Flores and Orozco, 2011).
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
#' @param fuz.num a non-negative \code{numeric}, the number of fuzzy
#' nucleosomes. Those nucleosomes are
#' distributed accordingly to selected distribution all over the region.
#' Those nucleosomes can overlap other well-positioned or fuzzy nucleosomes.
#'
#' @param fuz.var a non-negative \code{numeric}, the maximum variance of the
#' fuzzy nucleosomes. This variance
#' can be different than the one used for the well-positioned
#' nucleosome reads.
#'
#' @param max.cover a positive \code{numeric}, the maximum coverage for one
#' nucleosome. The final coverage can have a higher value than \code{max.cover}
#' since reads from differents nucleosomes can be overlapping. Default = 100.
#'
#' @param nuc.len a positive \code{integer}, the nucleosome length.
#' Default = 147.
#'
#' @param len.var a positive \code{numeric}, the variance of the distance
#' between a forward read and its paired reverse read. Default = 10.
#'
#' @param lin.len a non-negative \code{integer}, the length of the DNA linker
#' DNA. Default = 20.
#'
#' @param read.len a positive \code{integer}, the length of each of the
#' paired-end reads. Default = 40.
#'
#' @param rnd.seed a single value, interpreted as an \code{integer}, or
#' \code{NULL}. If a \code{integer} is given, the value is used to set the seed
#' of the random number generator. By fixing the seed, the generated results
#' can be reproduced. Default = \code{NULL}.
#'
#' @param distr the name of the distribution used to generate the nucleosome
#' map. The choices are : \code{"Uniform"}, \code{"Normal"} and
#' \code{"Student"}. Default = \code{"Uniform"}.
#'
#' @param offset a non-negative \code{integer}, the number of bases used to
#' offset all nucleosomes and reads. This is done to ensure that all
#' nucleosome positions and read alignment are of positive values.
#'
#' @return an \code{list} of \code{class} "syntheticSample" containing the
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
#' @author Pascal Belleau, Rawane Samb, Astrid Deschenes
#'
#' @examples
#'
#'## Generate a synthetic map with 20 well-positioned + 10 fuzzy nucleosomes
#'## using a Normal distribution with a variance of 30 for the well-positioned
#'## nucleosomes, a variance of 40 for the fuzzy nucleosomes and a seed of 15
#'## Because of the fixed seed, each time is going to be run, the results
#'## are going to be the seed
#'res <- syntheticNucReadsFromDist(wp.num = 20, wp.del = 0, wp.var = 30,
#'fuz.num = 10, fuz.var = 40, rnd.seed = 15, distr = "Normal",
#'offset = 1000)
#'
#' @export
syntheticNucReadsFromDist <- function(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                                    max.cover = 100,
                                    nuc.len = 147,
                                    len.var = 10,
                                    lin.len = 20,
                                    read.len = 40,
                                    rnd.seed = NULL,
                                    distr = c("Uniform", "Normal", "Student"),
                                    offset)
{
    ## Validation of parameters not used by syntheticNucMapFromDist() function
    syntheticNucReadsValidation(read.len, offset)

    ## Get call information
    cl <- match.call()

    ## Create a nucleosome map using the parameters of the user
    map <- syntheticNucMapFromDist(wp.num = wp.num, wp.del = wp.del,
                                    wp.var = wp.var, fuz.num = fuz.num,
                                    fuz.var = fuz.var, max.cover = max.cover,
                                    nuc.len = nuc.len, len.var = len.var,
                                    lin.len = lin.len, rnd.seed = rnd.seed,
                                    as.ratio = FALSE, distr = distr)

    return(createNucReadsFromNucMap(map, read.len, offset, cl))
}



#' @title Generate a synthetic nucleosome map containing forward and
#' reverse reads
#'
#' @description Generate a synthetic nucleosome map using a synthetic
#' nucleosome map.
#'
#' This function is using a modified version of the syntheticNucMap() function
#' from Bioconductor nucleR package (Flores and Orozco, 2011).
#'
#' @param syntheticNucMap a \code{list} of \code{class} "syntheticNucMap"
#'
#' @param read.len the length of each of the paired-end reads. Default = 40.
#'
#' @param offset the number of bases used to offset all nucleosomes and reads.
#' This is done to ensure that all nucleosome positions and read alignment
#' are of positive values.
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
#' @author Pascal Belleau, Rawane Samb, Astrid Deschenes
#'
#' @examples
#'
#'## Generate a synthetic map with 20 well-positioned + 10 fuzzy nucleosomes
#'## using a Normal distribution with a variance of 30 for the well-positioned
#'## nucleosomes, a variance of 40 for the fuzzy nucleosomes and a seed of 15
#'## Because of the fixed seed, each time is going to be run, the results
#'## are going to be the seed
#'
#'syntheticMap <- syntheticNucMapFromDist(wp.num = 20, wp.del = 0,
#'wp.var = 30, fuz.num = 10, fuz.var = 40,
#'rnd.seed = 335, as.ratio = FALSE, distr = "Uniform")
#'
#'res <- syntheticNucReadsFromMap(syntheticMap, read.len = 45, offset = 1000)
#'
#' @export
syntheticNucReadsFromMap <- function(syntheticNucMap, read.len = 40, offset) {

    ## Validate that object is of type "syntheticNucMap"
    if (!is(syntheticNucMap, "syntheticNucMap")) {
        stop("syntheticNucMap must be an object of class \'syntheticNucMap\'")
    }

    ## Validation of remaining parameters
    syntheticNucReadsValidation(read.len, offset)

    ## Get call information
    cl <- match.call()

    return(createNucReadsFromNucMap(syntheticNucMap, read.len, offset, cl))
}


#' @title Generate a graph of a synthetic nucleosome map
#'
#' @description Generate a graph for
#' a list marked as an \code{syntheticNucMap} class
#'
#' @param x a list marked as an \code{syntheticNucMap}  class
#'
#' @param ... \code{...} extra arguments passed to the \code{plot} function
#'
#' @examples
#'
#' ## Generate a synthetic map with 20 well-positioned nucleosomes, 5 fuzzy
#' ## nucleosomes and 10 deleted nucleosomes using a Student distribution
#' ## with a variance of 10 for the well-positioned nucleosomes,
#' ## a variance of 20 for the fuzzy nucleosomes and a seed of 15
#' syntheticMap <- syntheticNucMapFromDist(wp.num = 30, wp.del = 10,
#' wp.var = 10, fuz.num = 5, fuz.var = 20, rnd.seed = 15,
#' distr = "Student")
#'
#' ## Create graph using the synthetic map
#' plot(syntheticMap, xlab="Position", ylab="Coverage")
#'
#' @author Astrid Deschenes, Rawane Samb
#' @importFrom IRanges coverage
#' @importFrom graphics plot lines abline points legend polygon
#' @export
plot.syntheticNucMap <- function(x, ...) {

    ## Extract as.ratio parameter if present
    as.ratio = FALSE
    if (exists('syn.ratio', where=x)) {
        as.ratio = TRUE
    }

    ## Set Y axis maximum range
    max <- max(coverage(x$syn.reads), na.rm = TRUE) + 10

    ## Always set Y axis minimum to zero
    min <- 0

    ## Adapt values when as.ratio present
    if (as.ratio) {
        min <- min(x$syn.ratio, na.rm = TRUE)
    }

    # Plot coverage
    coverage <- c(0, as.integer(coverage(x$syn.reads)), 0)
    position <- c(0, 1:(length(coverage)-1))
    plot(position, coverage, type = "l", col = "gray",
         ylim = c(min, max), ...)
    polygon(c(0, position, 0), c(0, coverage, 0), col="gray", border = "gray")


    # Plot ratio, if asked for
    if (as.ratio) {
        lines(as.vector(x$syn.ratio), type = "l", col = "darkorange",
              lwd = 2)
        abline(h = 0, col = "black")
    }

    # Plot nucleosome positions
    space <- round(x$nuc.len/2)
    points(x$wp.starts  + space, x$wp.nreads,
            col = "forestgreen",  pch = 19)
    points(x$fuz.starts + space, x$fuz.nreads,
            col = "red", pch = 19)

    # Add legend
    if (as.ratio) {
        legend("top", c("Ratio", "Well", "Fuzzy", "Coverage"),
               fill = c("darkorange", "forestgreen", "red", "gray"),
               bty = "n", horiz = TRUE)
    } else {
        legend("top", c("Well", "Fuzzy", "Coverage"),
               fill = c("forestgreen", "red", "gray"), bty = "n",
               horiz = TRUE)
    }
}

#' @title Generate a graph of a synthetic nucleosome map using
#'
#' @description Generate a graph for
#' a list marked as an \code{syntheticNucReads} class
#'
#' @param x a list marked as an \code{syntheticNucReads} class
#'
#' @param ... \code{...} extra arguments passed to the \code{plot} function
#'
#' @examples
#'
#' ## Generate a synthetic map with 20 well-positioned nucleosomes, 5 fuzzy
#' ## nucleosomes and 10 deleted nucleosomes using a Student distribution
#' ## with a variance of 10 for the well-positioned nucleosomes,
#' ## a variance of 20 for the fuzzy nucleosomes and a seed of 15
#' syntheticNucSample <- syntheticNucReadsFromDist(wp.num = 30, wp.del = 10,
#' wp.var = 10, fuz.num = 5, fuz.var = 20, rnd.seed = 15,
#' distr = "Student", offset = 1000)
#'
#' ## Create graph using the synthetic map
#' plot(syntheticNucSample, xlab="Position", ylab="Coverage")
#'
#' @author Astrid Deschenes, Rawane Samb
#' @importFrom IRanges IRanges coverage
#' @importFrom graphics plot lines abline points legend polygon
#' @export
plot.syntheticNucReads <- function(x, ...) {

    ## Create a range using the sequences
    seqRanges <- IRanges(start=x$dataIP$start, end=x$dataIP$end)

    ## Set Y axis maximum range
    y_max <- max(coverage(seqRanges), na.rm = TRUE) + 10

    ## Always set Y axis minimum to zero
    y_min <- 0

    ## Set X axis minimum limit
    x_min <- min(x$dataIP$start) + 5

    ## Set X axis maximum limit
    x_max <- max(x$dataIP$end) + 5

    # Plot coverage
    coverage <- c(0, as.integer(coverage(seqRanges)), 0)
    position <- c(0, 1:(length(coverage)-1))
    plot(position, coverage, type = "l", col = "gray",
         ylim = c(y_min, y_max), xlim = c(x_min, x_max), ...)
    polygon(c(0, position, 0), c(0, coverage, 0), col="gray", border = "gray")


    # Plot nucleosome positions
    points(x$wp$nucleopos,  x$wp$nreads,  col = "forestgreen",  pch = 20)
    points(x$fuz$nucleopos, x$fuz$nreads, col = "red",          pch = 20)

    # Add legend
    legend("top", c("Coverage", "Well", "Fuzzy"),
               fill = c("gray", "forestgreen", "red"), bty = "n", horiz = TRUE)
}

