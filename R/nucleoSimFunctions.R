#' @title Generate a synthetic nucleosome map
#'
#' @description Generate a synthetic nucleosome map using the distribution
#' selected by the user. The user has choice between three different
#' distributions: Normal, Student and Uniform.
#'
#' This function is a modified version of the syntheticNucMap() function
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
#' nucleosomes can be overlapping. Default = 100.
#'
#' @param nuc.len a \code{numeric}, the nucleosome length. Default = 147.
#'
#' @param len.var a \code{numeric}, the variance of the distance between a
#' forward read and its paired reverse read. Default = 10.
#'
#' @param lin.len a \code{numeric}, the length of the DNA linker. Default = 20.
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
#' @param show.plot a \code{logical}, if \code{TRUE}, a plot of the coverage
#' map, with the nucleosome calls and optionally the calculated ratio is
#' outputed. Default = \code{FALSE}.
#'
#' @param distr the name of the distribution used to generate the nucleosome
#' map. The choices are : \code{"Uniform"}, \code{"Normal"} and
#' \code{"Student"}. Default = \code{"Uniform"}.
#'
#' @param ... arguments to be passed to \code{plot} function only
#' when \code{show.plot=TRUE}.
#'
#' @return an \code{list} of \code{class} "nucleoSim" containing the
#' following elements:
#' \itemize{
#'     \item \code{call} the matched call.
#'     \item \code{wp.starts} the start points of all well-positioned
#' nucleosomes.
#'     \item \code{wp.nreads} th number of repetitions of each well positioned
#' read.
#'     \item \code{wp.reads} a \code{IRanges} containing the well positioned
#' nucleosome reads.
#'     \item \code{fuz.starts} the start points of all the fuzzy nucleosomes.
#'     \item \code{fuz.nreads} the number of repetitions of each fuzzy
#' nucleosome read.
#'     \item \code{fuz.reads} a \code{IRanges} containing the fuzzy nucleosome
#' reads.
#'     \item \code{syn.reads} a \code{IRanges} containing all the synthetic
#' nucleosome reads (from fuzzy and well-positioned nucleosomes).
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
#'     fuz.num = 10, fuz.var = 40, show.plot = TRUE, rnd.seed = 15,
#'     distr = "Normal")
#'
#' ## Same output but without graph
#' syntheticNucMapFromDist(wp.num = 20, wp.del = 0, wp.var = 30,
#'     fuz.num = 10, fuz.var = 40, show.plot = FALSE,
#'     rnd.seed = 15, as.ratio = FALSE, distr = "Normal")
#'
#'
#' @author Rawane Samb, Astrid Deschenes
#' @importFrom stats runif rnorm rt
#' @importFrom IRanges IRanges
#' @importFrom graphics plot
#' @importFrom S4Vectors Rle
#' @importFrom IRanges coverage
#' @export
syntheticNucMapFromDist <- function(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                            max.cover = 100,
                            nuc.len = 147,
                            len.var = 10,
                            lin.len = 20,
                            rnd.seed = NULL,
                            as.ratio = FALSE,
                            show.plot = FALSE,
                            distr= c("Uniform", "Normal", "Student"), ...)
{
    # Get call information
    cl <- match.call()

    # Validate parameters
    syntheticNucMapFromDistValidation(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                            max.cover, nuc.len, len.var, lin.len, rnd.seed,
                            as.ratio, show.plot)

    # Validate distribution value
    distr <- match.arg(distr)

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

    # Delete some reads (set repetition time to 0)
    wp.nreads[round(runif(wp.del, min = 0, max = wp.num))] <- 0

    # Set each nucleosome as a repeated single start position
    wp.repstar <- rep(wp.starts, wp.nreads)

    # Add some variance to the starting points
    if (distr == "Normal")
    {
        wp.varstar <- wp.repstar + round(rnorm(length(wp.repstar), 0,
                                                wp.var^0.5))
        wp.varlen  <- nuc.len + round(rnorm(length(wp.repstar), 0,
                                                len.var^0.5))
    }
    else if (distr == "Student")
    {
        wp.varstar <- wp.repstar + round(rt(length(wp.repstar), 4))
        wp.varlen  <- nuc.len + round(rnorm(length(wp.repstar), 0,
                                           len.var^0.5))
    }
    else
    {
        wp.varstar <- wp.repstar + round(runif(length(wp.repstar),
                                                min=-wp.var,
                                                max=wp.var))
        wp.varlen <- nuc.len + round(rnorm(length(wp.repstar), 0,
                                           len.var^0.5))
    }

    # Putative reads
    #wp.reads <- IRanges(start=wp.varstar, width=nuc.len)
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
        fuz.varlen <- nuc.len + round(rnorm(length(fuz.repstar), 0,
                                                len.var^0.5))
    }
    else if (distr == "Student")
    {
        fuz.varstar <- fuz.repstar + round(rt(length(fuz.repstar), 4))
        fuz.varlen  <- nuc.len + round(rnorm(length(fuz.repstar), 0,
                                                len.var^0.5))
    }
    else
    {
        fuz.varstar <- fuz.repstar + round(runif(length(fuz.repstar),
                                                min = -fuz.var,
                                                max = fuz.var))
        fuz.varlen <- nuc.len + round(rnorm(length(fuz.repstar), 0,
                                                len.var^0.5))
    }
    # Overlapped reads
    fuz.reads <- IRanges(start = fuz.varstar, end = fuz.varstar + fuz.varlen)

    # ALL SYNTHETIC READS
    syn.reads <- c(wp.reads, fuz.reads)

    # RATIO AS HYBRIDIZATION (Tiling Array)
    if (as.ratio) {
        # Just put the same amount of reads as before randomly
        ctr.starts <- round(runif(length(syn.reads),
                                    min = 1,
                                    max = max(start(syn.reads))))

        # This time use a random read length, between 50 and 250
        ctr.widths <- round(runif(length(syn.reads), min = 50, max = 250))

        # "Control reads"
        ctr.reads <- IRanges(start = ctr.starts, width = ctr.widths)

        # ratio
        syn.ratio <- suppressWarnings(log2(as.vector(coverage(syn.reads))) -
                                        log2(as.vector(coverage(ctr.reads))))

        # Some lost bases... as reality
        syn.ratio[abs(syn.ratio) == Inf] <- NA
        syn.ratio <- Rle(syn.ratio)
    }

    # Preparing returned list
    result <- list(call = cl)

    result[["wp.starts"]] <- wp.starts
    result[["wp.nreads"]] <- wp.nreads
    result[["wp.reads"]]  <- wp.reads

    result[["fuz.starts"]] <- fuz.starts
    result[["fuz.nreads"]] <- fuz.nreads
    result[["fuz.reads"]]  <- fuz.reads

    result[["syn.reads"]]  <- syn.reads

    if(as.ratio) {
        result[["ctr.reads"]] <- ctr.reads
        result[["syn.ratio"]] <- syn.ratio
    }

    if(show.plot) {
        # Set Y-lim range
        max <- max(coverage(syn.reads), na.rm = TRUE)
        min <- 0
        if (as.ratio) {
            min <- min(syn.ratio, na.rm = TRUE)
        }

        # Plot main (coverage)
        plot(as.vector(coverage(syn.reads)), type = "h", col = "#AADDAA",
                            ylim = c(min, max), ...)

        # Plot ratio, if asked for
        if (as.ratio) {
            lines(as.vector(syn.ratio), type = "l", col = "darkorange",
                            lwd = 2)
            abline(h = 0, col = "darkorange4")
        }

        # Plot nucleosome positions
        points(wp.starts + 74, wp.nreads, col = "red", pch = 19)
        points(fuz.starts + 74, fuz.nreads, col = "blue", pch = 20)

        # Add legend
        if (as.ratio) {
            legend("top", c("Coverage", "Ratio", "Well-pos", "Fuzzy"),
                    fill = c("#AADDAA", "darkorange", "red", "blue"),
                    bty = "n", horiz = TRUE)
        } else {
            legend("top", c("Coverage", "Well-pos", "Fuzzy"),
                    fill = c("#AADDAA", "red", "blue"), bty = "n",
                    horiz = TRUE)
        }
    }

    return(result)
}



#' @title Generate a synthetic nucleosome map containing forward and
#' reverse reads
#'
#' @description Generate a synthetic nucleosome map using the distribution
#' selected by the user. The user has choice between three different
#' distributions: Normal, Student and Uniform.
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
#' nucleosomes can be overlapping. Default = 100.
#'
#' @param nuc.len the nucleosome length. Default = 147.
#'
#' @param len.var a \code{numeric}, the variance of the distance between a
#' forward read and its paired reverse read. Default = 10.
#'
#' @param lin.len the length of the DNA linker DNA. Default = 20.
#'
#' @param read.len the length of each of the paired reads. Default = 40.
#'
#' @param rnd.seed a single value, interpreted as an \code{integer}, or
#' \code{NULL}. If a \code{integer} is given, the value is used to set the seed
#' of the random number generator. By fixing the seed, the generated results
#' can be reproduced. Default = \code{NULL}.
#'
#' @param as.ratio a \code{logical}, if \code{TRUE},  a synthetic naked DNA
#' control map is created and the ratio between it and the nucleosome coverage
#' are calculated. It can be used to simulate hybridization ratio data, like
#' the one in Tiling Arrays. Both control map and calculated ratio are
#' returned. Default = \code{FALSE}.
#'
#' @param distr the name of the distribution used to generate the nucleosome
#' map. The choices are : \code{"Uniform"}, \code{"Normal"} and
#' \code{"Student"}. Default = \code{"Uniform"}.
#'
#' @param offset the number of bases used to offset all reads. This is done to
#' ensure that all nucleosome positions and read alignment are of positive
#' values.
#'
#' @return an \code{list} of \code{class} "nucleoSim" containing the
#' following elements:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{dataIP} a \code{data.frame} with all forward and reverse reads
#' for all well-positioned and fuzzy nucleosomes.
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
#'res <- syntheticSampleFromDist(wp.num = 20, wp.del = 0, wp.var = 30,
#'fuz.num = 10, fuz.var = 40, as.ratio = TRUE, rnd.seed = 15, distr = "Normal",
#'offset = 1000)
#'
#' @export
syntheticSampleFromDist <- function(wp.num, wp.del, wp.var, fuz.num, fuz.var,
                                max.cover = 100,
                                nuc.len = 147,
                                len.var = 10,
                                lin.len = 20,
                                read.len = 40,
                                rnd.seed = NULL,
                                as.ratio = FALSE,
                                distr = c("Uniform", "Normal", "Student"),
                                offset)
{
    ## Validate that offset is a non-negative integer
    if (!isInteger(offset) || offset < 0) {
        stop("offset must be a non-negative integer")
    }

    ## Get call information
    cl <- match.call()

    ## Create a nucleosome map using the parameters of the user
    map <- syntheticNucMapFromDist(wp.num = wp.num, wp.del = wp.del,
                                    wp.var = wp.var, fuz.num = fuz.num,
                                    fuz.var = fuz.var, max.cover = max.cover,
                                    nuc.len = nuc.len, len.var = len.var,
                                    lin.len = lin.len, rnd.seed = rnd.seed,
                                    as.ratio = as.ratio, show.plot = FALSE,
                                    distr = distr)

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
    ## Forward and reverse reads are generated forme read paired information
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

    # Create dataframe with well-positioned nucleosomes shifted by offset value
    wp <- data.frame(map$wp.starts + offset + round(nuc.len/2), map$wp.nreads)
    colnames(wp) <- c("nucleopos", "nreads")

    # Create dataframe with fuzzy nucleosomes shifted by offset value
    fuz <- data.frame(map$fuz.starts + offset + round(nuc.len/2),
                            map$fuz.nreads)
    colnames(fuz) <- c("nucleopos", "nreads")

    # Create returned result list
    result <- list(call = cl, dataIP = dataIP, wp = wp, fuz = fuz,
                   paired = paired)

    return(result)
}
