#' @title Gramian matrix, norms, and inner products of splines
#' @param Sp \code{Splinets} object;
#' @param norm_only logical, indicates if only the square norm of the 
#' elements in the input object is calculated; The default is \code{norm_only=FALSE};
#' @param sID vector of integers, the indicies specifying splines in the \code{Splinets} 
#' list \code{Sp} to be evaluated; If \code{sID=NULL} (default), then the inner products for all the pairs taken from the object are evaluated.
#' @param Sp2 \code{Splinets} object, the optional second  \code{Splinets}-object; The inner products between 
#' splines in \code{Sp} and in \code{Sp2} are evaluated, i.e. the cross-gramian matrix.
#' @param s2ID vector of integers, the indicies specifying splines in the \code{Sp2} to be considered in the cross-gramian;
#' @return 
#' \itemize{
#' \item \code{norm_only=FALSE} -- the Gram matrix of inner products of the splines within the input \code{Splinets}-objects is returned, 
#' \item \code{Sp2 = NULL} -- the non-negative definite matrix of the inner products of splines in \code{Sp} is returned,
#' \item both \code{Sp} and \code{Sp2} are non-\code{NULL} and contain splines \eqn{S_i}'s and \eqn{T_j}'s, respectively == 
#' the cross-gramian matris of the inner products for the pairs of splines \eqn{(S_i,T_j)} is returned,
#' \item \code{norm_only=FALSE}-- the vector of the norms of \code{Sp} is returned.
#' }
#' @description The function performs evaluation of the matrix of the inner products 
#' \eqn{\int S(t) \cdot T(t) dt }{\int S(t) * T(t) dt} of all the pairs of splines \eqn{S}, \eqn{T} from the input object.
#' The program utilizes the Taylor expansion of splines, see the reference for details.
#' @details If there is only one input \code{Splinet}-object, then the non-negative symmetrix matrix of the splines in this object is returned. 
#' If there are two input \code{Splinet}-objects, then the \eqn{m \times r}{m x r} matrix of the cross-inner product is returned, where \eqn{m} is 
#' the number of splines in the first object and \eqn{r} is their number in the second one. 
#' If only the norms are evaluated (\code{norm_only= TRUE}) it is always evaluating the norms of the first object. 
#' In the case of two input \code{Splinets}-objects, they should be over the same set of knots and of the same smoothness order. 
#' @export
#' @inheritSection Splinets-class References
#' @seealso \code{\link{lincomb}} for evaluation of a linear combination of splines;
#'  \code{\link{project}} for projections to the spaces of Splines; 
#' @example R/Examples/ExGramian.R
#' 

gramian = function(Sp , norm_only = FALSE, sID = NULL, Sp2 = NULL ,s2ID = NULL){
  knots = Sp@knots
  k = Sp@smorder
  S = Sp@der
  supp = Sp@supp
  
  n = length(knots)-2 #the number of internal knots
  d = length(S)       #the number of splines in the object
  full_supp = matrix(c(1,n+2), ncol = 2) # 1 x 2 matrix for the full support case
  
  if(length(supp) == 0) supp = rep(list(full_supp), d) #adding the list of the full supports corresponding to the list of splines in the Splinet object
  if(is.null(sID)) sID = 1:d
  for(i in sID) S[[i]] = sym2one(S[[i]], supp = supp[[i]])
  
  n_so = length(sID)
  FF = FF_Matrix(knots, k)
  D = as.matrix(1/1:(2*k+1))
  
  if(norm_only){
    Gram_M = numeric(length(sID)) #The vector of the norms will be returned
    for(i in 1:length(sID)){
      Gram_M[i] = inner_engine(S[[sID[i]]], S[[sID[i]]],
                               supp[[sID[i]]], supp[[sID[i]]], FF, D, n, k)
    }
  }else{#actual gramian matrices 
    if(is.null(Sp2)){#the single Splinets input
      Gram_M = matrix(numeric(n_so^2), n_so)
      for(i in 1:n_so){
        for(j in i:n_so){ #starting from i to avoid double computation in the single Splinets input case
          Gram_M[i,j] = inner_engine(S[[sID[i]]], S[[sID[j]]],
                                     supp[[sID[i]]], supp[[sID[j]]], FF, D, n, k)
        }
      }
      Gram_M = t(Gram_M-diag(diag(Gram_M),nrow=n_so)) + Gram_M #to obtain the second half below the diagonal in the single Splinets input case
    }else{#The case of two Splinets
      if(prod(knots != Sp2@knots) | (k != Sp2@smorder)){stop("The order or the knots of the two input splinets do not agree.")}
      S2 = Sp2@der
      supp2 = Sp2@supp
      d2 = length(S2)       #the number of splines in the object
      full_supp2 = matrix(c(1,n+2), ncol = 2) # 1 x 2 matrix for the full support case
      
      if(length(supp2) == 0) supp2 = rep(list(full_supp2), d2) #adding the list of the full supports corresponding to the list of splines in the Splinet object
      if(is.null(s2ID)) s2ID = 1:d2
      for(i in s2ID) S2[[i]] = sym2one(S2[[i]], supp = supp2[[i]]) #trasnsforming the derivative matrices to the one sided representation
      n_so2 = length(s2ID)
      Gram_M = matrix(0,nrow=n_so, ncol= n_so2)
      for(i in 1:n_so){
        for(j in 1:n_so2){ 
          Gram_M[i,j] = inner_engine(S[[sID[i]]], S2[[s2ID[j]]],
                                     supp[[sID[i]]], supp2[[s2ID[j]]], FF, D, n, k)
        }
      }
    }
    }#The end of main calculation
    
  return(Gram_M)
}