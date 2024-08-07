#' @title Evaluating splines at given arguments.
#' @description For a \code{Splinets}-object \code{S} and a vector of arguments \code{t}, 
#' the function returns the matrix of values for the splines in \code{S}. The evaluations are done 
#' through the Taylor expansions, so on the \eqn{i}th interval for  
#' \eqn{t\in [\xi_i,\xi_{i+1}]}{[\xi_i,\xi_[i+1]]}:
#' \deqn{S(t)=\sum_{j=0}^{k} s_{i j} \frac{(t-\xi_{i})^j}{j!}.}{S(t) = \sum(l=1:k) (t-\xi_i)^l * 1/l! * s_il.}
#' For the zero order splines which are discontinuous at the knots, the following convention is taken. 
#' At the LHS knots the value is taken as the RHS-limit,  at the RHS knots as the LHS-limit. 
#' The value at the central knot for the zero order and an odd number of knots case is assumed to be zero. 
#' @param object \code{Splinets} object;
#' @param sID vector of integers, the indicies specifying splines in the \code{Splinets} list 
#' to be evaluated; If \code{sID=NULL}, then all splines in the \code{Splinet}-object are evaluated. The default value
#' is \code{NULL}.  
#' @param x vector, the arguments at which the splines are evaluated; If \code{x} is
#' \code{NULL}, then the splines are evaluated over regular grids per each interval of the support. The default value is \code{x=NULL}.
#' @param N integer, the number of points per an interval between two consequitive knots at which the splines are evaluated.
#' The default value is \code{N = 250};
#'
#' @return The \code{length(x) x length(sID+1)} matrix containing the argument values, in the first column, 
#' then, columnwise, values of the subsequent splines.
#' @export
#' 
#' @inheritSection Splinets-class References
#'
#' @seealso \code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects;
#' \code{\link{plot,Splinets-method}} for plotting \code{Splinets}-objects; 
#' @example R/Examples/ExEvspline.R
#' 
#' 
evspline = function(object, sID = NULL, x = NULL, N = 250){
  xi = object@knots
  n = length(xi) - 2
  supp = object@supp
  S = object@der
  l = length(S) #the number of splines in the object
  
  if(length(supp) == 0){ #the case when the whole range is the support. 
    supp = rep(list(matrix(c(1,n+2), ncol = 2)), l)
  }
  
  if(is.null(sID)) sID = 1:l #all splines from the input
  
  r = range(xi) 
  if(is.null(x)){
    if(object@equid){
      x <- seq(r[1], r[2], length.out = N*(n+2)) #regular grid
      }else{
      x <- NULL
      for(i in 1:(n+1)) x = c(x, seq(xi[i], xi[i+1], length.out = N))
      }
    }else{
    rr = range(x)
    if(rr[2] > r[2] | rr[1] < r[1]){
      stop("The range of 'x' is not covered by the full range of knots")
    }
  }
  
  y = matrix(0, nrow = length(x), ncol = length(sID))
  
  for(i in 1:length(sID)){
    y[,i] = evaluate_spline(xi, supp[[sID[i]]], object@smorder, S[[sID[i]]], x)
  }
  evsp=cbind(x,y)
  return(unname(evsp))
}
