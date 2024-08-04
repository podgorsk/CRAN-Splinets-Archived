#' @title B-splines, periodic B-splines and their orthogonalization
#' 
#' @description The B-splines (periodic B-splines)  are either given in the input or generated inside the routine. Then, given 
#' the B-splines and the argument \code{type}, the routine additionally generates a \code{Splinets}-object
#' representing an orthonormal spline basis obtained from a certain 
#' orthonormalization of the B-splines. Orthonormal spline bases are obtained by one of the following methods:
#' the Gram-Schmidt method, the two-sided method, and/or the splinet algorithm, which is the default method.
#' All spline bases are kept in the format of \code{Splinets}-objects.
#' @param knots \code{n+2} vector, the knots (presented in the increasing order); It is not needed, when
#' \code{Bsplines} argumment is not \code{NULL}, in which the case the knots from \code{Bsplines} are inherited.
#' @param smorder integer,  the order of the splines, the default is \code{smorder=3}; Again it is inherited from the
#' \code{Bsplines} argumment if the latter is not \code{NULL}.
#' @param type string, the type of the basis; The following choices are available 
#' \itemize{
#'   \item \code{'bs'} for the unorthogonalized B-splines,
#'   \item \code{'spnt'} for the orthogonal splinet (the default),
#'   \item \code{'gsob'} for the Gramm-Schmidt (one-sided) O-splines,
#'   \item \code{'twob'} for the two-sided O-splines.
#'  } 
#' @param Bsplines \code{Splinet}-object, the basis of the B-splines (if not \code{NULL}); 
#' When this argument is not \code{NULL} the first two arguments  
#' are not needed since they will be inherited from \code{Bsplines}.
#' @param norm logical, a flag to indicate if the output B-splines should be normalized;
#' @param periodic logical, a flag to indicate if B-splines will be of periodic type or not;
#' @return Either a list \code{list("bs"=Bsplines)} made of a single \code{Splinet}-object \code{Bsplines} 
#' when \code{type=='bs'}, which represents the B-splines (the B-splines are normalized or not, depending
#' on the \code{norm}-flag), or a list of two \code{Splinets}-objects: \code{list("bs"=Bsplines,"os"=Splinet)}, 
#' where \code{Bsplines} are either computed (in the input \code{Bspline= NULL}) or taken from the input \code{Bspline}
#' (this output will be normalized or not depending on the \code{norm}-flag),
#' \code{Splinet} is the B-spline orthognalization determined by the input argument \code{type}. 
#' @details 
#'  The B-spline basis, if not given in 
#' the input, is computed 
#' from  the following recurrent (with respect to the smoothness order of the B-splines) formula
#' \deqn{
#' B_{l,k}^{\boldsymbol \xi }(x)=
#' \frac{x- {\xi_{l}}
#'  }{
#' {\xi_{l+k}}-{\xi_{l}}
#' }
#' B_{l,k-1}^{\boldsymbol \xi}(x)
#' +
#'  \frac{{\xi_{l+1+k}}-x }{ {\xi_{l+1+k}}-{\xi_{l+1}}}
#'  B_{l+1,k-1}^{\boldsymbol \xi}(x), l=0,\dots, n-k.
#' }{
#'  B_lk(x)=(x-\xi_l)/(\xi_{l+k}-\xi_l) * B_{lk-1}(x)
#'  +
#'  (\xi_{l+1+k}-x)/(\xi_{l+1+k}-\xi_{l+1}) * B_{l+1k-1}(x), l=0,\dots, n-k
#'  } 
#'  The dyadic algorithm that is implemented takes into account efficiencies due to the equally space knots 
#' (exhibited in the Toeplitz form of the Gram matrix) only if the problem is fully dyadic, i.e. if the number of 
#' the internal knots is \code{smorder*2^N-1}, for some integer \code{N}. To utilize this efficiency it may be advantageous, 
#' for a large number of equally spaced knots, to choose them so that their number follows the fully dyadic form.
#' An additional advantage of the dyadic form is the complete symmetry at all levels of the support. The algorithm works with 
#' both zero boundary splines and periodic splines. 
#' @export
#' @inheritSection Splinets-class References
#' @example R/Examples/ExSplinet.R 
#' @seealso \code{\link{project}} for projecting into the functional spaces spanned by the spline bases; 
#' \code{\link{lincomb}} for evaluation of a linear combination of splines;
#' \code{\link{seq2dyad}} for building the dyadic structure for a splinet of a given smoothness order;
#' \code{\link{plot,Splinets-method}} for visualisation of splinets; 


splinet = function(knots=NULL, smorder = 3, type = 'spnt', Bsplines=NULL, periodic= FALSE,  norm=F){
  if(!is.null( Bsplines)){
   periodic= Bsplines@periodic
  }
  if (!periodic) {
    splnt=splinet1(knots=knots, smorder = smorder, type = type, Bsplines=Bsplines,  norm=norm)  
  }else{
    splnt=splinet2(knots=knots, smorder = smorder, type = type, Bsplines=Bsplines,  norm=norm) 
  }
  return(splnt)
}
