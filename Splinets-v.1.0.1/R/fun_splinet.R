#' @title B-splines and their orthogonalization
#' 
#' @description The B-splines are either given in the input or generated inside the routine. Then, given 
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
#' An additional advantage of the dyadic form is the complete symmetry at all levels of the support. 
#' @export
#' @inheritSection Splinets-class References
#' @example R/Examples/ExSplinet.R 
#' @seealso \code{\link{project}} for projecting into the functional spaces spanned by the spline bases; 
#' \code{\link{lincomb}} for evaluation of a linear combination of splines;
#' \code{\link{seq2dyad}} for building the dyadic structure for a splinet of a given smoothness order;
#' \code{\link{plot,Splinets-method}} for visualisation of splinets; 


splinet = function(knots=NULL, smorder = 3, type = 'spnt', Bsplines=NULL, norm=FALSE){
  
  #------------------------------#
  # S1: generating bspline basis #
  #------------------------------#
  if(!is.null(Bsplines)){ #inheriting the arguments if B-splines are in the input
    knots=Bsplines@knots
    smorder=k=Bsplines@smorder
    n = length(knots) - 2
    
  }else{
    k = smorder
    n = length(knots) - 2
    
    #In the case knots are not in the increasing order they are sorted
    if(min(diff(knots))<0){
      knots=sort(unique(knots))
      cat("Knots were not given in the strictly increasing order, which is required.\n
          Ordered  knots with removed ties are replacing the input values.\n")
    }
    #
    
    #Creating a generic 'Splinets' object to store computed spline bases
    so = new("Splinets", knots = knots, smorder = k) #it checks among other things if knots are equidistant
    #and sets 'so@equid' to a proper value, see `setClass` in 
    # Splinets-class definition
    so@type = "bs"
    supp = cbind(1:(n-k+1),(1:(n-k+1))+k+1)
    supp = so@supp = lapply(seq_len(nrow(supp)), function(i) supp[i,,drop=FALSE]) #seq_len(n)=1:n, 
    #creates a list of matrices not vectors
    if(so@equid){
      xi = knots[1:(k+2)] #so there will be only computations to get one element of the basis
    } else{               #the remaining elements of the basis will have identical derivative matrices
      xi = knots
    }
    n = length(xi)-2
    
    S = list()
    S[[1]] = array(rep(1, n+1), dim = c((n+1), 1))
    for(ord in 1:k){
      S[[ord+1]] = array(0, dim = c((ord+1)*(n-ord+1), ord+1))
      l1 = dim(S[[ord]])[1]; l2 = dim(S[[ord]])[2]
      TS = array(rep(0, (l1+l1/ord)*l2), dim = c(l1+l1/ord, l2)) # temp augmented S matrix
      TS[-(ord+1)*(1:(n-ord+2)), ] = S[[ord]]
      TS1 = head(TS, -(ord+1)); TS2 = head(tail(TS, -ord),-1)
      c1 = coeff1(xi, ord); c2 = coeff2(xi, ord)
      lam1 = lambda1(xi, ord); lam2 = lambda2(xi, ord)
      S[[ord+1]][,1] = c1*lam1*TS1[,1] + c2*lam2*TS2[,1]
      S[[ord+1]][,ord+1] = c1*ord*TS1[,ord] + c2*ord*TS2[,ord]
      if(ord>1){
        for(j in 2:ord){
          S[[ord+1]][,j] = c1*((j-1)*TS1[,j-1] + lam1*TS1[,j]) + c2*((j-1)*TS2[,j-1] + lam2*TS2[,j])
        }
      }
    }
    S = S[[k+1]]
    # Transform to the symmetric representation of derivative matrices
    SS = list()
    n = length(knots)-2
    if(so@equid){
      for(i in 1:(n-k+1)){
        SS[[i]] = sym2one(rbind(S, numeric(k+1)), supp[[i]], inv = TRUE)
      }
    } else{
      for(i in 1:(n-k+1)){
        SS[[i]] = sym2one(rbind(S[(i+(i-1)*k):(i+i*k),], numeric(k+1)), supp[[i]], inv = TRUE)
      }
    }
    so@der = SS
    n_so = length(SS)
    Bsplines=so  #The Bsplines are computed 
    } 
  #=======================
  #The end of the B-splines 
  #------------------------------#
  # S2: normalization of bspline basis #
  #------------------------------#  
  # normalization of the B-splines (one does not assume that the input splines are orthogonalized)
  a = array(1/sqrt(gramian(Bsplines, norm_only = TRUE)), dim=c((n-k+1),1))
  so = Bsplines
  for(i in 1:length(a)){
    so@der[[i]] = a[i]*so@der[[i]]
  }
  if(norm==TRUE){Bsplines=so} #Normalization of the output B-splines
  splnt=list("bs"=Bsplines) #The B-spline part of the output list
  n_so=length(so@der) #The number of osplines in the basis
  #------------------------------#
  # S3: Orthogonalization of B-splines #
  #------------------------------#  
  if(type == 'gsob'){
    H = bandmatrix(knots, k, so@der, so@supp)
    P = grscho(diag(n_so), H) #the generic algorithm for GS orthogonalization
    so = lincomb(so, t(P))
    so@type = 'gsob'
    splnt$os=so
  }
  # S2.2: twosided obasis
  if(type == 'twob'){
    H = bandmatrix(knots, k, so@der, so@supp)
    P = sgrscho(diag(n_so), H) #the generic algorithm for symmetric GS orthogonalization
    so = lincomb(so, t(P))
    so@type = 'twob'
    splnt$os=so
  }
  # S2.3: splinet
  if(type == 'spnt'){
    
    H=bandmatrix(knots,k,so@der,so@supp)
    
    P=dyadiag(H,k,so@equid) #The main band-matrix diagonalization algorithm
    so = lincomb(so, t(P))
    n=dim(H)[1] + k - 1
    N=log2((n+1)/k)
    if(N-floor(N)!=0){so@type = "spnt"}else{so@type = "dspnt"} #flagging if not dyadic
    
    splnt$os=so 
  }
  
  return(splnt)
}

