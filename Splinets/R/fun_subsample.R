#' @title Subsampling from a set of splines
#' @description The function constructs a \code{Splinets}-object that is made of subsampled
#' elements of the input \code{Splinets}-object.
#' The input objects have to be of the same order and over the same knots.
#' @param Sp \code{Splinets}-object, a collection of \code{s} splines;
#' @param ss vector of integers, the coordinates from \code{1:s};
#' @return An \code{Splinets}-object containing \code{length(ss)} splines that are selected from the input object./
#' @export
#' @details The output \code{Splinet}-object made of subsampled splines is always is of the regular type, i.e. SLOT \code{type='sp'}.
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects;
#' \code{\link{construct}} for constructing such a \code{Splinets}-object;
#' \code{\link{gather}}  for combining \code{Splinets}-objects;
#' \code{\link{refine}} for refinment of a spline to a larger number of knots;
#' \code{\link{plot,Splinets-method}} for plotting \code{Splinets}-objects;
#' @example R/Examples/ExSubsample.R
#'
#'
subsample=function(Sp,ss){
  ss=as.vector(ss) #to handle the case of a single number to be treated as a vector
  SPout=Sp
  n1=length(ss)

  SPout@der=list()
  SPout@supp=list()
  if(length(Sp@supp)!=0){
    for(i in 1:n1){
      SPout@supp[[i]]=Sp@supp[[ss[i]]] #Subsampling the supports 
    }
  }
  for(i in 1:n1){
    SPout@der[[i]]=Sp@der[[ss[i]]] #Subsampling the derivative matrix over the support
  }
  SPout@type='sp'
  return(SPout)

} #The end of the function


