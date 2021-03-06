% Generated by roxygen2 (4.0.2): do not edit by hand
\name{multsims}
\alias{multsims}
\title{Share simulations}
\usage{
multsims(nsims, names, nmons, reps, ndays, PCs, keeps, cms, sds,
  unequal = NULL, days = NULL, cut = 1, thres = pi/4, prnt = F,
  sderr = 0.01)
}
\arguments{
\item{nsims}{number of iterations}

\item{names}{vector of source names (e.g. c("traffic", "fireworks", "soil")) corresponding to PCs}

\item{nmons}{number of monitors}

\item{reps}{number of monitors per subregion}

\item{ndays}{number of observations}

\item{PCs}{positive part of PC from sample data}

\item{keeps}{share info for creating data}

\item{cms}{vector of lognormal means for sources}

\item{sds}{vector of lognormal sds for sources}

\item{unequal}{vector of numbers to switch subregions}

\item{days}{vector of days for each monitor}

\item{cut}{cutoff for eigenvalues (see nmsource), default is 1.}

\item{thres}{threshold for share angle cutoff}

\item{prnt}{Print simulation iteration (default = F)}
}
\description{
\code{multsims} Performs share simulation for multiple datasets
}
\details{
This is a function to compare results from SHARE and mAPCA for
multiple simulated datasets
}

