#' @title Linear transformation of splines.
#'
#' @description A linear combination of the splines \eqn{S_j} in the input object is computed according to  
#' \deqn{R_i=\sum_{j=0}^{d} a_{i j} S_j,\, i=1,\dots, l.}
#' and returned as a \code{Splinet}-object.
#' @param object \code{Splinets} object containing \code{d} splines;
#' @param A \code{l x d} matrix; coefficients of the linear transformation,
#' @param reduced logical; If \code{TRUE} (default), then the linear combination is 
#' calculated accounting for the actual support sets (recommended for sparse splines), 
#' if \code{FALSE}, then the 
#' full support computations are used (can be faster for lower dimension or non-sparse cases).
#' @param SuppExtr logical; If \code{TRUE} (default), the true support is extracted, otherwise, full range 
#' is reported as the support. Applies only to the case when \code{reduced=FALSE}.   
#'
#' @return A \code{Splinet}-object that contains \code{l} splines obtained by linear combinations of  
#' using coefficients in rows of \code{A}. The  \code{SLOT type} of the output splinet objects is \code{sp}.
#' 
#' 
#' @export
#' @inheritSection Splinets-class References
#' @seealso \code{\link{exsupp}} for extracting the correct support; 
#' \code{\link{construct}} for building a valid spline; 
#' \code{\link{rspline}} for random generation of splines;
#' @example R/Examples/ExLincomb.R
#' 
#' 

lincomb = function(object, A, reduced = TRUE, SuppExtr = TRUE){
  k = object@smorder
  S = object@der
  xi = object@knots
  supp = object@supp
  n = length(xi)-2
  if(dim(A)[2] != length(S)){
    stop("The number of columns of input A matrix is different from the number of splines in the input object")
  }
  
  if(reduced){
    # do linear combination using reduced der matrix
    temp = lincomb_supp(S, supp, A, n)
  } else{
    # do linear combination using full der matrix
    temp = lincomb_full(S, supp, k, n, A, SuppExtr)
  }
  
  object@type = "sp"
  object@der = temp$S
  object@supp = temp$supp
  
  return(object)
}