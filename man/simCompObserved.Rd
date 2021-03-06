% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branchCorr_source.r
\name{simCompObserved}
\alias{simCompObserved}
\title{Simulate discretely observed data from a stochastic compartmental model}
\usage{
simCompObserved(t.end, initPopulation, rates, progStructure, obsTimes,
  vecOutput = FALSE, maxEvents = 999999999999)
}
\arguments{
\item{t.end}{The total time until end of simulation}

\item{initPopulation}{A vector of length equal to total number of compartments/types; entries contain the initial population of each type}

\item{rates}{A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
(renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)}

\item{progStructure}{A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends}

\item{obsTimes}{A vector containing the observation times}

\item{vecOutput}{Logical, whether to return a vectorized representation of the matrix of observation times by population counts}

\item{maxEvents}{The maximum number of events to simulate}
}
\value{
A matrix of population sizes of each compartment at each observaiton time. Rows index cell type; columns index observation times. Returns
the vectorized form of thie matrix if vecOutput=TRUE.
}
\description{
This code simulates data analogously to \code{\link{simCompartments}}, but records populations at a specified list of observation times.
}
\examples{
progStructure <- c(1,1,1,2,2)
initPopulation <- c(10,0,0,0,0,0,0,0)
obsTimes <- 1:5
rates <- c(.3,.2,.5,.06, .03, 2, 4, 5, 3, 1, .15,.5,.8,.1,.05)
simCompObserved(5, initPopulation, rates,progStructure,obsTimes)
}

