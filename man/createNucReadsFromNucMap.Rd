% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nucleoSimInternFunctions.R
\name{createNucReadsFromNucMap}
\alias{createNucReadsFromNucMap}
\title{Create a synthetic nucleosome reads from a synthetic nucleosome map}
\usage{
createNucReadsFromNucMap(map, read.len, offset, call)
}
\arguments{
\item{map}{a \code{list} of \code{class} "syntheticNucMap"}

\item{read.len}{the length of each of the paired-end reads. Default = 40.}

\item{offset}{the number of bases used to offset all nucleosomes and reads.
This is done to ensure that all nucleosome positions and read alignment
are of positive values.}

\item{call}{the function call.}
}
\value{
an \code{list} of \code{class} "syntheticNucReads" containing the
following elements:
\itemize{
\item \code{call} the matched call.
\item \code{dataIP} a \code{data.frame} with the chromosome name, the
starting and ending positions and the direction of all forward
and reverse reads for all well-positioned and fuzzy nucleosomes.
\item \code{wp} a \code{data.frame} with the positions of all the
well-positioned nucleosomes, as well as the number of paired-reads
associated to each one.
\item \code{fuz} a \code{data.frame} with the positions of all the fuzzy
nucleosomes, as well as the number of paired-reads associated to each one.
\item \code{paired} a \code{data.frame} with the starting and ending
positions of the reads used to generate the paired-end reads.
}
}
\description{
Generate a synthetic nucleosome map using a synthetic
nucleosome map.
}
\examples{

## Generate a synthetic map with 30 well-positioned nucleosomes, 5 fuzzy
## nucleosomes and 6 deleted nucleosomes using a Student distribution
## with a variance of 10 for the well-positioned nucleosomes,
## a variance of 15 for the fuzzy nucleosomes and a seed of 1335
map_call <- call("syntheticNucMapFromDist", wp.num = 30, wp.del = 6,
wp.var = 10, fuz.num = 5, fuz.var = 15, rnd.seed = 1335,
distr = "Student")
syntheticMap <- eval(map_call)

syntheticReads <- nucleoSim:::createNucReadsFromNucMap(syntheticMap,
read.len = 40, offset = 1000, call = map_call)

}
\author{
Astrid Deschenes
}
\keyword{internal}
