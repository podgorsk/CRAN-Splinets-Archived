#' @title Derivatives of splines
#' @description The function generates a \code{Splinets}-object which contains the first order
#' derivatives of all the splines from the input \code{Splinets}-object.
#' The function also verifies the support set of the output to provide the accurate information about 
#' the support sets by excluding regions over which the original function is constant. 
#' @param object \code{Splinets} object of the smoothness order \code{k};
#' @param epsilon positive number, controls removal of knots from the support; If the derivative is smaller than this number, it is considered 
#' to be zero and the corresponding knots are removed from the support.The default value is \code{1e-7}.  
#' @return A \code{Splinets}-object of the order \code{k-1} that also contains the updated information about the support set.
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{integra}} for generating the indefinite integral of a spline that can be viewed 
#' as the inverse operation to \code{deriva}; 
#' \code{\link{dintegra}} for the definite integral of a spline;
#' @example R/Examples/ExDeriva.R
#' 
#'  
deriva = function(object, epsilon = 1e-7){
  k = object@smorder
  
  S = object@der
  supp = object@supp
  xi = object@knots
  
  n = length(xi)-2
  d = length(S) #The number of splines in the object. 
  
  full_supp = matrix(c(1,n+2), ncol = 2)
  if(length(supp) == 0){#It means that the full support is assumed for each spline in the object
    supp = rep(list(full_supp), d)
  }
  
  if(k==0){stop("The zero order case, the derivative is trivially equal to zero.")}

  
  for(i in 1:d){
    S[[i]] = S[[i]][,-1,drop=FALSE] #removing the first column is equivalent to taking the derivative
    
      if(sum(apply(S[[i]] == 0, 1, prod)) > 0){ # fix constant part, i.e. if there some rows in the matrix
                                              # of derivatives that are zero it requires a correction of support
                                              # this correspond to the original function to be constant over some stretch.
      temp = exsupp(S[[i]], supp[[i]]) #extracting (correcting) the support
      supp[[i]] = temp$rsupp 
      S[[i]] = temp$rS
    }
  }
  
  object@der = S
  object@smorder = k-1
  object@taylor = object@taylor[, -(k+1), drop = FALSE]
  object@supp = supp
  object@type = 'sp' #derivative is always of 'sp' type by default
  return(object)
}