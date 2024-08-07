#' @title Organizing indices in a spline basis in the net form
#' 
#' @description This auxiliary function generates the map between the sequential order and the dyadic net structure
#' of a spline basis. It works only with indices so it can be utilized to any basis in the space of 
#' splines with the zero-boundary conditions. The function is useful for creating the dyadic structure of the graphs
#' and whenever a reference to the \code{k}-tuples and the levels of support is needed. 
#' @param n_sp positive integer, the number of splines to be organized into the dyadic net; 
#' The dyadic net does not need to be fully dyadic, i.e. \code{n_sp} does not need to be equal to \eqn{k2^n-1}, 
#' where \eqn{n} is the number of the internal knots. See the references for more details.
#' @param k the size of a tuple in the dyadic net; It naturally corresponds to the smoothness order of splines for which the 
#' net is build. 
#' @return The double indexed list of single row matrices of positive integers in the range \code{1:n_sp}.
#' Each vector has typically the length \code{k} and some of them may correspond to incomplete tuplets and thus can be 
#' shorter. The first index in the list points to the level in the dyadic structure, the second one to the 
#' the number of the tuplet at the given level. The integers in the vector pointed by the list 
#' correspond to the sequential index of the element belonging to this tuplet. 
#' @inheritSection Splinets-class References
#' @export
#' @seealso \code{\link{plot,Splinets-method}} for plotting splinets in the dydadic graphical representation; 
#' \code{\link{lincomb}} for evaluation of a linear combination of splines;
#' \code{\link{refine}} for refinment of a spline to a larger number of knots; 
#' @example R/Examples/ExSeq2dyad.R 
#' 
seq2dyad = function(n_sp, k){
  res=net_structure(n_sp,k)
  colnames(res) <- NULL
  net=list()
  for(i in 1:max(res[,2]))
  {
    net[[i]]=list()
    ind=res[,2]==i #the indexes of splines that belong to level i
    ind_tupl=sort(unique(res[ind,3])) #the indexes of tuplets in the level i
    for(j in 1:length(ind_tupl)){
      aa=res[,3]==ind_tupl[j] #flags the indices that belong to the j-th tuplet at the level i
      net[[i]][[j]]=matrix(sort(res[aa,1]),nrow=1)
    }
  }
  return(net)
}
