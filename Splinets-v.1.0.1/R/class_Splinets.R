#' @title The class to represent a collection of splines
#' @description The main class in the \code{splinets}-package used for representing a collection of splines.
#' @slot knots numeric \code{n+2} vector, a vector of n+2 knot locations presented in the increasing order and without ties;
#' @slot smorder non-negative integer, the smoothnes order of the splines, i.e. the highest order of non-zero derivative;
#' @slot equid logical, indicates if the knots are equidistant;
#' Some computations in the equidistant case are simpler so this information helps to account for it.
#' @slot supp list (of matrices), \itemize{
#'   \item \code{length(supp)==0} -- the full support set for all splines,
#'   \item \code{length(supp)==N} -- support sets for \code{N} splines;
#' } If non-empty, a list containing \code{Nsupp x 2} matrices (of positive integers). 
#' If \code{Nsupp} is equal to one it should be a row matrix (not a vector).
#' The rows in the matrices, \code{supp[[i]][l,]}, \code{l in 1:Nsupp} represents the indices of the knots that are the endpoints of the intervals in the support sets. 
#' Each of the support set is represented as a union of disjoint \code{Nsupp} intervals, with knots as the endpoints. Outside the set (support), the spline vanishes. 
#' Each matrix in this list is ordered so the rows closer to the top correspond to the intervals closer to the LHS end of
#' the support. 
#' @slot der list (of matrices); a list of the length \code{N} containing
#' \code{sum(supp[[i]][,2]-supp[[i]][,1]+1) x (smorder+1)} matrices, where \code{i} is the index running through the list.
#' Each matrix in the list includes the values of the derivatives at the knots in the support of the corresponding spline.
#' @slot taylor \code{(n+1) x (smorder+1)}, if \code{equid=FALSE}, or \code{1 x (smorder+1)}
#' if \code{equid=TRUE},  columnwise vectors  of the Taylor expansion coefficients at the knots; 
#' Vectors instead of matrices are recognized properly.
#' The knot and order dependent matrix of rows of coefficients used in the Taylor expansion of splines.
#' Once evaluated it can be used in computations for any spline of the given order over the given knots.
#' The columns of this matrix are used for evaluation of the values of the splines in-between knots,
#' see the references for further details.
#' @slot type string, one of the following character strings: \code{bs},\code{gsob},\code{twob},\code{dspnt},\code{spnt},\code{sp}; The default is \code{sp} which 
#' indicates any unstructured collection of splines. The rest of the strings indicate different \emph{spline bases}: 
#' \itemize{
#'   \item \code{bs} for B-splines,
#'   \item \code{gsob} for Gram-Schmidt O-splines,
#'   \item \code{twob} for two-sided O-splines,
#'   \item \code{dspnt} for a fully dyadic splinet,
#'   \item \code{spnt} for a non-dyadic splinet.
#' }  
#' @slot epsilon numeric (positive), an accuracy used to detect a problem with the conditions
#' required for the matrix of the derivatives (controls relative deviation from the conditions);
#'
#' @return running \code{new("Splinets")} return an object that belongs to the class \code{Splinets}, with the initialization of the default
#' values for the fields.
#' @export
#' @section References:
#' Liu, X., Nassar, H., Podg\eqn{\mbox{\'o}}{o}rski, K. (2019) "Splinets -- efficient orthonormalization of the B-splines."  <arXiv:1910.07341>.
#'
#' Podg\eqn{\mbox{\'o}}{o}rski, K. (2021) 
#' "\code{Splinets} -- splines through the Taylor expansion, their support sets and orthogonal bases." <arXiv:2102.00733>.
#' 
#' @seealso \code{\link{is.splinets}} for evaluation of a \code{Splinets}-object; \code{\link{construct}} for constructing a \code{Splinets}-object; 
#' \code{\link{plot,Splinets-method}} for plotting methods for \code{Splinets}-objects;
#'
#' @example R/Examples/ExSplinetsObject.R
#' @importFrom methods callNextMethod new

setClass("Splinets",
         representation(
           knots="vector", smorder="numeric", equid="logical",
           supp="list", der="list", taylor = "matrix",type = "character", epsilon="numeric"
         ),
         prototype(
           knots=c(0,1), smorder=0, equid=FALSE, der=list(as.matrix(c(1,1))), type="sp", epsilon=0.0000001
           )
         #prototype sets the initial values for slots, if desired 
)

setMethod("show","Splinets", print.class) #This will use function 'print.class' to show an object from the class

setMethod("initialize", "Splinets", function(.Object, ...) { #This will be used for the function 'new()'
  .Object <- callNextMethod() #It is not obvious what this is supposed to do but it is a standard, it somehow calls 
                              #the currently described method on '.Object' after finishing this method definition.
                              #Not a very precise description but it works. 
  n <- length(.Object@knots) - 2
  k <- .Object@smorder
  
  # 1) validate knots
  if(n < k){
    stop(paste("SLOT 'knots' in a 'splinets' object should be a vector of at least length ",k+2,"\n",
               "Reconsider a vector of increasing values for SLOT 'knots'."))
  }
  m <- min(diff(.Object@knots))
  if(m<=0){
    .Object@knots=sort(unique(.Object@knots))
    cat("NOTE: Knots are not in the strictly increasing order, which is required.\n
        The knots have been ordered and ties have been removed.")
    n=length(.Object@knots)-2
  }
  
  # 2) set 'equid' attr. 
  h=diff(.Object@knots)
  errkn=max(h)-min(h)
  .Object@equid <- (errkn < 0.01*.Object@epsilon) #if the distance between knots is approximately the same
  
  # 3) set 'taylor'
  if(min(dim(as.matrix(.Object@taylor)))==1){
    if(k!=0){#one dimensional vector has to be treated as a column matrix
      .Object@taylor=matrix(.Object@taylor,nrow=1)
    }else{ #one dimensional vector has to be treated as a row matrix
      .Object@taylor=matrix(.Object@taylor,ncol=1)
    }
  }#To assure that illogical treatment of vectors and matices in R is fixed.
  
  if(dim(.Object@taylor)[2]!=(k+1)){ #If the matrix has been evaluated before, then it has 'k+1' columns
    if(.Object@equid==TRUE){
      .Object@taylor<-taylor_coeff(.Object@knots[1:2],k) #all the coefficients for all other knots are the same as for
      #the first two
    }else{
      .Object@taylor<-taylor_coeff(.Object@knots,k)
    }
  }
  
  .Object
} )


