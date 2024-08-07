#' @title Indefinite integrals of splines
#' 
#' @description The function generates the indefinite integrals for given input splines. 
#' The integral is a function 
#' of the upper limit of the definite integral and is a spline of the higher order that does not satisfy the zero boundary conditions at the RHS endpoint,
#' unless the definite integral over the whole range is equal to zero. 
#' Moreover, the support of the function is extended in the RHS up to the RHS end point 
#' unless the definite integral of the input is zero, in which the case the support is extracted from the obtained spline. 
#' @param object a \code{Splinets} object of the smoothness order \code{k};
#' @param epsilon non-negative number indicating accuracy when close to zero value are detected; This accuracy is used in 
#' when the boundary conditions of the integral are checked.  
#' @return A \code{Splinets}-object with order \code{k+1} that contains the indefinite integrals
#' of the input object.
#' @details  
#' The value on the RHS is not zero, so the zero boundary condition typically is not satisfied and the support is 
#' is extended to the RHS end of the whole domain of splines. However, the function returns proper support
#' if the original spline is a derivative of a spline that satisfies the boundary conditons. 
#' @export
#' @inheritSection Splinets-class References
#' @seealso \code{\link{deriva}} for computing derivatives of splines; \code{\link{dintegra}} for the definite integral;
#' @example R/Examples/ExIntegra.R
#' 
#' 
integra = function(object , epsilon=1e-07){
  xi = object@knots
  k = object@smorder
  supp = object@supp
  S = object@der
  newS = S #the list to keep the new entries of the matrices of the derivatives for integral

  taylor = object@taylor
  d = length(S)
  n = length(xi)-2
  
  if(object@equid){
    taylor = matrix(rep(taylor, n+1), byrow = TRUE, nrow = n+1)
    taylor = cbind(taylor, diff(xi)^(k+1)/factorial(k+1) )
  } else{
    taylor = cbind(taylor, diff(xi)^(k+1)/factorial(k+1) ) #computing SLOT taylor for the integral
  }
  
  
  if(length(supp) == 0){ #To assign the support intervals to the case of the full supports
    full_supp = matrix(c(1,n+2), ncol = 2)
    supp = rep(list(full_supp), d)
  }
  
  #Throughout the next loop the partial support can be extracted in the case that the integral at the last 
  #knot is zero. If the support at all splines is the full support, then the output SLOT should be 
  #an empty list to control it the following indicator is set
  
  FullSupp=TRUE
  
  for(i in 1:d){ #running through the splines in the object
    temps=supp[[i]] #support components for the ith spline
    Nsupp=dim(temps)[1] #the number of components
    
    res = int_engine(supp = temps, taylor = taylor, S = S[[i]]) #computing the derivative matrices over each component
    
    
    S[[i]] = sym2one(res,temps) #At this point the derivatives do not account for the constant value between the 
                    #constant values between support components neither between the last component and the end of 
                    #the knots range if there are any knots after the last component. 
  
    
    nS = matrix(0,nrow=(n+2),ncol=(k+2)) #The space for the derivative values that assumes the full range support
    
    
    compS=der_split(S[[i]],temps) #This allows below to aivod counting indices by separating the matrices corresponding 
                                 #to different support components (the components will be in one-sided representation)
    
    if(Nsupp>1){#if there is more than one component in the support
      for(j in 1:(Nsupp-1)){#running through the original support components
      
        nS[temps[j,1]:temps[j,2],]=compS[[j]] #the j-th component is added
        #First the constant value between the components of the support
        #in between of the original support ends constant value is assigned
        nS[(temps[j,2]+1):(temps[j+1,1]-1),1]=nS[temps[j,2],1] #it will repeat the value from the last knot
     }
    }  
       nS[temps[Nsupp,1]:temps[Nsupp,2],]=compS[[Nsupp]] #the last component is added
       #The values at the knots after the last component
       if(temps[Nsupp,2]<(n+2)){
        nS[(temps[Nsupp,2]+1):(n+2),1]=nS[temps[Nsupp,2],1]
       } #This puts constant value to the RHS of the last support component.
  
      #Next to extract the support components if the value at the end is zero
      if(abs(nS[n+2,1])<epsilon){#the case of the spline (the definite integral of the object is zero)
        
        #In the next the symmetric form of the matrix must be given 
        
        aa=exsupp(sym2one(nS,inv=TRUE)) #extracting support from the full support matrix
        if(aa$rsupp[1,1]!=1 | aa$rsupp[1,2]!=(n+2)){#the returned support is not the full support
          
          FullSupp=FALSE #output SLOT 'supp' will not be the empty list (if at least one of the splines 
                         #has smaller than the full support)
          supp[[i]]=aa$rsupp #new partial support
          nS=sym2one(aa$rS,aa$rsupp) #corresponding matrix of the derivatives, onesided at this point
        
        }
        
    }#The end of extracting support
    if(FullSupp==TRUE){
      newS[[i]]=sym2one(nS,inv=TRUE) #
      }else{
      newS[[i]]=sym2one(nS,supp[[i]],inv=TRUE) #one needs to get back to the symmetric representation
      }
  }
  
  

  #Expanding to the full support since typically the integral will have full support 
  if(FullSupp==TRUE){supp=list()} #This is how the full support is marked in the object
  
  object@supp=supp
  object@der = newS
  object@smorder = k+1
  

  
  object@type = 'sp' #integral is always of 'sp' type by default
  
  if(object@equid){
    object@taylor = taylor[1,,drop=FALSE]
  }else{
    object@taylor = taylor
  }
  
  return(object)
}