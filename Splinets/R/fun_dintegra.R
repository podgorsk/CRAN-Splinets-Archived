#' @title Calculating the definite integral of a spline. 
#' 
#' @description The function calculates the definite integrals of the splines in an input  \code{Splinets}-object.
#' @param object \code{Splinets}-object;
#' @param sID vector of integers, the indicies specifying for which splines in the \code{Splinets}-object 
#' the definite integral is to be evaluated; If \code{sID=NULL}, then the definite integral of all splines in the object 
#' are calculated. The default is \code{NULL}. 
#'
#' @return A \code{length(sID) x 2} matrix, with the first column holding the id of splines and the second 
#' column holding the corresponding definite integrals. 
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{integra}} for generating the indefinite integral; 
#' \code{\link{deriva}} for generating derivative functions of splines;
#' @example R/Examples/ExDintegra.R
dintegra = function(object, sID = NULL){
  xi = object@knots
  S = object@der
  supp = object@supp
  taylor = object@taylor
  
  if(object@equid){
    taylor = matrix(rep(taylor, n+1), byrow = TRUE, nrow = n+1)
  }
  
  n = length(xi)-2
  d = length(S)
  full_supp = matrix(c(1,n+2), ncol = 2)
  if(length(supp) == 0){
    supp = rep(list(full_supp), d)
  }
  if(is.null(sID)){
    sID = 1:d
  }
  res = cbind(sID, numeric(length(sID)))
  colnames(res) = c("Spline ID", "dIntegral")
  for(i in 1:length(sID)){
    res[i, 2] = dint_engine(knots = xi, supp = supp[[sID[i]]], taylor = taylor, S = S[[sID[i]]])
  }
  return(res)
}