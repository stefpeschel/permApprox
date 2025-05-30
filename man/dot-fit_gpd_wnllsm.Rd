% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dot_gpd_fit_methods.R
\name{.fit_gpd_wnllsm}
\alias{.fit_gpd_wnllsm}
\title{WNLSM and WNLLSM parameter estimation for the GPD}
\usage{
.fit_gpd_wnllsm(
  x,
  boundary = NULL,
  eval_point = NULL,
  method = "WNLLSM",
  tol = 1e-08,
  shapeIni = 0.01,
  scaleIni = 0.01,
  F0 = 0,
  plotopt = FALSE,
  ...
)
}
\arguments{
\item{x}{data vector}

\item{boundary}{numeric. Maximum value at which the GPD probability density
function must be positive. If \code{NULL}, the maximum of the data vector
\code{x} is used.}

\item{tol}{numeric giving the desired accuracy of the optimization
process.}

\item{shapeIni, scaleIni}{initial values for the shape and scale parameters,
i.e., values where the optimization starts.}
}
\description{
Weighted nonlinear least squares moments (WNLSM) estimation and
  likelihood WNLSM (WNLLSM) estimation for the two-parameter Generalized
  Pareto Distribution (GPD) proposed by \cite{Zhao et al. (2019)}.
}
\details{
A simplified version of the R code has been provided by
  \cite{Zhao et al. (2019)}. The function has been adapted to ensure that
  the GPD density is non-zero at boundary + eps or max(x) + eps.
}
\references{
\insertRef{Zhao2019new}{permApprox}
}
\keyword{internal}
