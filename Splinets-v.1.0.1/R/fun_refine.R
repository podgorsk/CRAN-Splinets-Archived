#' @title Refining splines through adding knots
#' 
#' @description Any spline of a given order remains a spline of the same order if one considers it on a bigger set of knots than the original one.
#'              However, this embedding changes the \code{Splinets} representation of the so-refined spline. 
#'              The function evaluates the corresponding \code{Splinets}-object. 
#' @param object \code{Splinets}-object, the object to be represented as a \code{Splinets}-object over a refined set of knots;
#' @param mult positive integer, refining rate; The number of the knots to be put equally spaced between the existing knots.
#' @param newknots \code{m} vector, new knots; 
#' The knots do not need to be ordered and knots from the input \code{Splinets}-object knots are allowed since any ties are resolved. 
#' @return A \code{Splinet} object with the new refined knots and the new matrix of derivatives is evaluated at the new knots combined with the original ones. 
#' @details The function merges new knots with the ones from the input \code{object}. It utilizes \code{deriva()}-function to evaluate the derivative at the refined knots. 
#' It removes duplications of the refined knots, and account also for the not-fully supported case.
#' In the case when the range of the additional knots extends beyond the knots of the input \code{Splinets}-object,
#' the support sets of the output \code{Splinets}-object account for the smaller than the full support. 
#' @export
#' @seealso \code{\link{deriva}} for computing derivatives at selected points; 
#' \code{\link{project}} for an orthogonal projection into a space of splines;
#' 
#' @inheritSection Splinets-class References
#' @example R/Examples/ExRefine.R
#' 

refine = function(object,  mult=2, newknots=NULL){
k = object@smorder
S = object@der
xi = object@knots
supp = object@supp
n = length(xi)-2
d = length(S) #The number of splines in the object. 

newobject=object           #We start with the input object which will be subsequently modified into the output
newobject@equid=FALSE      #It may become TRUE if the original object is such and newknots=NULL, see below
if(is.null(newknots)){#Thus new knots are governed by 'mult'
  if(object@equid==TRUE){ #equidistant case
    newobject@equid=TRUE
    newobject@knots=seq(xi[1],xi[n+2],by=(xi[2]-xi[1])/mult)
    newxi=newobject@knots
    newobject@taylor=taylor_coeff(newxi[1:2],k) #all the coefficients for all other knots are the same as for
    #the first two
  }else{
      newxi=seq(xi[1],xi[2]-(xi[2]-xi[1])/mult,by=(xi[2]-xi[1])/mult)
      for(j in 1:n){newxi=c(newxi,seq(xi[j+1],xi[j+2]-(xi[j+2]-xi[j+1])/mult,by=(xi[j+2]-xi[j+1])/mult))}
      newxi=c(newxi,xi[n+2])
      newobject@knots=newxi
      newobject@taylor=taylor_coeff(newxi,k)
    }
}else{ #it is disregarding 'mult' parameter
  newobject@knots=unique(sort(c(xi,newknots)))
  newobject@taylor=taylor_coeff(newobject@knots,k)
}


newxi=newobject@knots         

newobject@type = "sp" #when we refine the resulting object always becomes of the type 'sp', i.e. an unspecified collection of splines



newn=length(newxi)-2

derlist=list() #the list of derivatives
derlist[[1]]=object
for(i in 1:k) derlist[[i+1]]=deriva(derlist[[i]]) #computing the derivatives of of the input Splines up to the k-th order

if(length(supp) == 0 & prod(range(newxi)==range(xi))){#It means that the full support is assumed for each spline in the object
                                                      #The case when the range of the refinement is the same as the original range of knots
   newS=matrix(0,ncol=k+1,nrow=newn+2)
   for(r in 1:d){#running through the splines in the object
     for(i in 1:(k+1)){#running through the columns of the matrix of derivatives (one-sided version for the highest derivative)
       drv=derlist[[i]]
       evsp=evspline(drv,sID = r ,x=newxi)
       newS[,i]=evsp[,2,drop=FALSE] 
     }
     newS=sym2one(newS,inv=TRUE) #the way the evaluation is performed on the the zero order splines makes it one-sided values at the last column
     newobject@der[[r]]=newS
     newobject@supp=list() #the full support case
   }
}else{
  if(length(supp) == 0){#the case of the smaller range of the input than the refined splines
    for(i in 1:d){
      supp[[i]]=matrix(c(1,n+2),ncol=2) #the full support explicitly because the output will not have the full support
    }
  }
  newobject@supp=supp #The new object will have the same size of the support but indexes will be changed below
  for(i in 1:d){#running through the splines in the object
    #first, the new support indexes has to be provided based on the previous ones
    ns=dim(supp[[i]])[1] #the number of the support components
    for(j in 1:ns){
      newobject@supp[[i]][j,1]=which(newxi==xi[supp[[i]][j,1]]) #assigning corrected indexes for the support - LEFT
      newobject@supp[[i]][j,2]=which(newxi==xi[supp[[i]][j,2]]) #RIGHT
    }
    newobject@der[[i]]=matrix(0,nrow=sum(newobject@supp[[i]][,2]-newobject@supp[[i]][,1]+1),ncol=k+1)
  } #the new supports has been created and empty matrid for derivatives of proper size
  #Evaluation of the matrices of derivatives
  for(i in 1:d){ # running through the splines in the object
    lastend=0
    for(j in 1:dim(newobject@supp[[i]])[1]){#running through the components of the support
      newS=matrix(0,ncol=k+1,nrow=newobject@supp[[i]][j,2]-newobject@supp[[i]][j,1]+1) #creating a matrix for the derivatives 
      locknots=newxi[newobject@supp[[i]][j,1]:newobject@supp[[i]][j,2]] #local new knots to have the derivatives evaluated
      for(l in 1:(k+1)){#running through the columns of the matrix of derivatives (one-sided version for the highest derivative)
        drv=derlist[[l]]
        evsp=evspline(drv,sID = i,x=locknots)
        newS[,l]=evsp[,2,drop=FALSE] 
      }#the end of evaluations of the derivatives on a given component
      newS=sym2one(newS,inv=TRUE) #the way the evaluation is performed on the the zero order splines makes it one-sided values at the last column
      curbeg=lastend+1
      curend=lastend+newobject@supp[[i]][j,2]-newobject@supp[[i]][j,1]+1
      newobject@der[[i]][curbeg:curend,]=newS #assigning the j-th component
      lastend=curend #this many entries in the matrix has been filled out. 
     }#the evalutions for the components in the support set
  }#the evaluations for the splines in the object
}
return(newobject)
}
