#' @title Combining two \code{Splinets} objects
#' @description The function returns the \code{Splinets}-object that gathers two input \code{Splinets}-objects together.
#' The input objects have to be of the same order and over the same knots.
#' @param Sp1 \code{Splinets} object;
#' @param Sp2 \code{Splinets} object;
#' @return \code{Splinets} object, contains grouped splines from the input objects;
#' @export
#' 
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{is.splinets}} for diagnostic of the \code{Splinets}-objects;
#' \code{\link{construct}} for constructing such an object;
#' \code{\link{subsample}}  for subsampling \code{Splinets}-objects;
#' \code{\link{plot,Splinets-method}} for plotting \code{Splinets}-objects;
#' @example R/Examples/ExGather.R
#'
#'
gather=function(Sp1,Sp2){
   SP=Sp1
   if(length(SP@supp)!=0){
      if(length(Sp2@supp)!=0){
         SP@supp=c(SP@supp,Sp2@supp)
      }else{
         v=list(t(as.matrix(c(1,length(Sp2@knots)))))
         sz=length(Sp2@der)
         SP@supp=c(SP@supp,rep(v,sz))     
      }
   }else{
      if(length(Sp2@supp)!=0){
         v=list(t(as.matrix(c(1,length(SP@knots)))))
         sz=length(Sp1@der)
         SP@supp=c(rep(v,sz),Sp2@supp)   
      }
   }
   SP@der=c(Sp1@der,Sp2@der)
   return(SP)
   
} #The end of the function