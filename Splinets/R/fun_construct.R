#' @title Construction of a \code{Splinets} object
#'
#' @description The function constructs a \code{Splinets} object correspond to a single spline (size=1)
#' from a vector of knots and a matrix of proposed derivatives.
#' The matrix is tested for its correctness like in \code{is.splinets} and adjusted using one of the implemented methods.
#' @param knots \code{n+2} vector, the knots over which the spline is built; 
#' There should be at least \code{2*smorder+4} of knots.
#' @param smorder integer, the order of smoothness;
#' @param matder \code{(n+2)x(smorder+1)} matrix, the matrix of derivatives; 
#' This matrix will be corrected if does not correspond to a proper spline.
#' @param supp vector, either empty or two integers representing the single interval support; 
#' @param mthd string, one of the three methods for correction of the matrix of derivative:
#' \describe{
#'  \item{\code{'CRLC'}}{matching mostly the highest derivative,}
#'  \item{\code{'CRFC'}}{matching mostly the function values at the knots,}
#'  \item{\code{'RRM'}}{balanced matching between all derivatives;}
#' }
#' The default method is \code{'RRM'}, see the paper on the package for further details about the methods.
#' @return A \code{Splinets}-object corresponding to a single spline.
#' @details The function constructs a \code{Splinet}-object only over a single interval support. 
#' Combining with the function \code{lincom} allows to introduce a multi-component support.  
#' @export
#' 
#' @inheritSection Splinets-class References
#'
#' @seealso \code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects;
#' \code{\link{gather}} and \code{\link{subsample}}  for combining and subsampling \code{Splinets}-objects, respectively,
#' \code{\link{plot,Splinets-method}} for a plotting method for \code{Splinets}-objects;
#' \code{\link{lincomb}} for combining splines with more complex than a single interval support sets; 
#' @example R/Examples/ExConstruct.R
#'
#'
construct=function(knots,smorder,matder,supp=vector(),mthd='RRM'){
   
   if(length(supp)==0){#The full support is marked by an empty list
      support=list()
      }else{
     if(min(dim(t(supp)))==1){#we consider only one component supports in this construction
        support=list(matrix(supp,ncol=2)) #for the list of the matrices form of the support
     }else{
        stop("The input 'supp' is not representing a one-component support, which is required.")
        } 
     if(dim(matder)[1]!=(support[[1]][1,2]-support[[1]][1,1]+1)){
        stop('The size of the matrix of the derivative does not match the support size.')
        }
   }
   
   der=list()
   der[[1]]=matder
   cat("\nUsing  method",mthd,"to correct the derivative matrix entries.\n")
   build=new("Splinets",knots=knots,smorder=smorder,supp=support,der=der)
   build=is.splinets(build)
   build=build[[2]]
   Vd=build@taylor
   if(mthd!='RRM'){
     if(length(support)!=0){Vd=Vd[build@supp[[1]][1,1]:(build@supp[[1]][1,2]-1),]}#These are Taylor's coefficient between the knots of the support
     build@der[[1]]=correct(matder,Vd,mthd)
   }
   cat("\nThe matrix derivative is now corrected by method",mthd,".\n") 
    return(build)

  } #The end of the function


