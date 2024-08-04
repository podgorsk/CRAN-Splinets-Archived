#' @title Data on tire responses to a rough road profile
#'
#' @description These are simulated data of tire responses to a rough road at the high-transient
#' event. The simulations have been made based on the fit of the so-called Slepian model
#' to a non-Gaussian rough road profile. Further details can be found in the reference. The
#' responses provided are measured at the wheel and thus describing the tire response. 
#' There are 100 functional measurments, kept column-wise in the matrix.
#' Additionally, the time instants of the measurements are given as the first column in the matrix.
#' Since the package uses the so-called "lazy load", the matrix 
#' is directly available without an explicit load of the data.
#' This means that \code{data(tire)} does not need to be invoked.
#' The data were saved using \code{compress='xz'} option, which requires 3.5 or higher version of R.
#' The data are uploaded as a dataframe, thus \code{as.matrix(tire)} is needed if the matrix form is required. 
#' @docType data
#'
#' @usage data(tire)
#'
#' @format numerical \code{4095 x 101} dataframe: \code{tire} 
#'
#' @keywords datasets
#'
#' @seealso \code{\link{truck}} for a related dataset;
#' @section References:
#' Podg\eqn{\mbox{\'o}}{o}rski, K, Rychlik, I. and Wallin, J. (2015)
#'  Slepian noise approach for gaussian and Laplace moving average processes. 
#'  Extremes, 18(4):665–695, <doi:10.1007/s10687-015-0227-z>. 
#'
#' @example R/Examples/ExTruckTire.R
#' 
"tire"



