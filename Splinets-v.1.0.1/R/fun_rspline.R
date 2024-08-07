#' @title Random splines
#' @description The function simulates a random \code{Splinets}-object that is made of random splines with the center 
#' at the input spline and the matrix of derivatives has the added error term of the form 
#' \deqn{
#' \boldsymbol \Sigma^{1/2}\mathbf Z \boldsymbol \Theta^{1/2},
#'  } 
#' where  \eqn{\mathbf Z} is a \eqn{(n+2)\times (k+1)} matrix having iid standard normal variables
#' as its entries, while \eqn{\boldsymbol \Sigma} and \eqn{\boldsymbol \Theta} are matrix parameters. 
#' This matrix error term is then corrected by one of the methods and thus resulting in a matrix of derivatives at knots corresponding to a valid spline. 
#' @param S \code{Splinets}-object with \code{n+2} knots and of the order of smoothness \code{k}, representing the center of randomly simulated splines; 
#' When the number of splines in the object is bigger than one, only the first spline in the object is used. 
#' @param N positive integer, size of the sample;
#' @param Sigma matrix; \itemize{
#' \item If \code{(n+2)x(n+2)} matrix, it controls correlations between derivatives of the same order at different knots.
#' \item If a positive number, it represents a diagonal \code{(n+2)x(n+2)} matrix with this number on the diagonal.
#' \item If a \code{n+2} vector, it 
#' represents a diagonal \code{(n+2)x(n+2)} matrix with the vector entries on the diagonal.
#' \item If \code{NULL} (default) represents the identity matrix.
#' } 
#' @param Theta matrix; 
#' \itemize{
#' \item If  \code{(k+1)x(k+1)}, this controls correlations between different derivatives at each knot.
#' \item If a positive number, it represents a diagonal matrix with this number on the diagonal.
#' \item If a \code{k+1} vector, it 
#' represents a diagonal matrix with the vector entries on the diagonal.
#' \item If  \code{NULL} (default), it represents the \code{k+1} identity matrix;
#' }
#' @param mthd string, one of the three methods: RCC, CR-LC, CR-FC, to adjust random error matrix so it corresponds to a valid spline;
#'
#' @return A \code{Splinets}-object that contains \code{N} generated splines constituting an iid sample of splines;
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{is.splinets}} for diagnostics of the \code{Splinets}-objects;
#' \code{\link{construct}} for constructing a \code{Splinets}-object;
#' \code{\link{gather}}  for combining \code{Splinets}-objects into a bigger object;
#' \code{\link{subsample}}  for subsampling \code{Splinets}-objects;
#' \code{\link{plot,Splinets-method}} for plotting \code{Splinets}-objects;
#' @example R/Examples/ExRspline.R
#' @importFrom stats convolve rnorm

rspline=function(S,N=1,Sigma=NULL,Theta=NULL,mthd='RRM'){

invisible(capture.output(chk<-is.splinets(S))) #checking if the mean value is a Splinets object
                                     #suppressing the message (chk=) would not work must be (chk<-)
if(chk[[1]]==FALSE){stop("The mean value parameter is not a Splinets object.")}

  n=length(S@knots)[1]-2

  k=dim(S@der[[1]])[2]-1
if(is.null(Sigma)){sqS=diag(n+2) #standard independent rows
}else{
  Sigma=as.matrix(Sigma) 
  mS=min(dim(Sigma)[1],dim(Sigma)[2])
  MS=max(dim(Sigma)[1],dim(Sigma)[2])
  
  if(MS==1){
    sqS=diag(sqrt(Sigma[1,1]),n+2) #iid rows case
  }
  if((MS!=mS)&&mS==1){sqS=diag(as.vector(sqrt(Sigma)))} #independent rows case
  if(MS==mS&& mS>1){
    spS=eigen(Sigma,symmetric=TRUE)
    sqS=spS$vectors%*%diag(sqrt(spS$values))%*%t(spS$vectors) #square roots of covariances
  }
}
if(is.null(Theta)){
  sqT=diag(k+1) #standard independent columns
}else{
  Theta=as.matrix(Theta)
  mT=min(dim(Theta)[1],dim(Theta)[2])
  MT=max(dim(Theta)[1],dim(Theta)[2])
  if(MT==1){sqT=diag(sqrt(Theta[1,1]),k+1)} #iid columns case
  if((MT!=mT)&&mT==1){sqT=diag(as.vector(sqrt(Theta)))} #independent rows case
  if(MT==mT&& mT>1){
    spT=eigen(Theta,symmetric=TRUE)
    sqT=spT$vectors%*%diag(sqrt(spT$values))%*%t(spT$vectors) #square roots of covariance
  }
}
for(i in 1:N)
{
  Z=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
  derT=sqS%*%Z%*%sqT
  if(length(S@supp)==0){
    derT=derT+S@der[[1]]
  }else{
    B=1 #Counting cumulatively knots in the support
    for(j in 1:dim(S@supp[[1]])[1]){ #running through support intervals
      E=S@supp[[1]][j,2]-S@supp[[1]][j,1]+B
      derT[S@supp[[1]][j,1]:S@supp[[1]][j,2],]=derT[S@supp[[1]][j,1]:S@supp[[1]][j,2],]+S@der[[1]][B:E,]
      B=E+1
    }
  }

invisible(capture.output(rspl<-construct(S@knots,k,derT,mthd=mthd)))
if(i>1){rspls=gather(rspls,rspl)}else{rspls=rspl}
}
return(rspls)
} #The end of the function

