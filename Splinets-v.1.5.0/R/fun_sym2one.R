#' @title Switching between representations of the matrices of derivatives
#'
#' @description A technical but useful transformation of the matrix of derivatives form the one-sided 
#' to symmetric representations, or a reverse one. It allows for switching between the standard representation of the matrix 
#' of the derivatives for \code{Splinets} which is symmetric around the central knot(s) to the one-sided that yields 
#' the RHS limits at the knots, which is more convenient for computations. 
#' @param S \code{(m+2) x (k+1)} numeric matrix, the derivatives in one of the two representations; 
#' @param supp \code{(Nsupp x 2)} or \code{NULL} matrix, row-wise the endpoint indices of the support intervals; If it 
#' is equal to \code{NULL} (which is also the default), then the full support is assumed. 
#' @param inv logical; If \code{FALSE} (default), then the function assumes that the input is 
#' in the symmetric format and transforms it to the left-to-right format. If \code{TRUE}, then 
#' the inverse transformation is applied.
#' @details The transformation essentially changes only the last column in \code{S}, i.e. the highest (discontinuous) derivatives so that 
#' the one-sided representation yields the right-hand-side limit. 
#' It is expected that the number of rows in \code{S} is the same as the total size of the support
#' as indicated by \code{supp}, i.e. if \code{supp!=NULL}, then \code{sum(supp[,2]-supp[,1]+1)=m+2}. 
#' If the latter is true, than all derivative submatrices of the components in \code{S} will be reversed.
#' However, this condition formally is not checked in the code, which may lead to switch of 
#' the representations only for parts of the matrix \code{S}.
#' @return A matrix that is the respective transformation of the input.
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{Splinets-class}} for the description of the \code{Splinets}-class; 
#' \code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects; 
#' @example R/Examples/ExSym2one.R
#'
#'
sym2one = function(S, supp=NULL, inv=FALSE){
  h = dim(supp)[1]
  if(h == 1 || is.null(h) ){ #The single interval support (including the full support)
    S = sym2one_single(S, inv = inv) 
  }else{
    dd = supp[,2]-supp[,1] # for extracting der matrix given support
    st = 1 # start row of der
    for(j in 1:h){
      en = st + dd[j]
      S[st:en, ] = sym2one_single(S[st:en, ,drop=FALSE], inv = inv)
      st = en+1
    }
  }
  return(S)
}
