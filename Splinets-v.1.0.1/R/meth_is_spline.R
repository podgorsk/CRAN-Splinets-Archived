#' @title Diagnostics of splines and their generic correction
#'
#' @description The method performs verification of the properties of SLOTS of an object belonging to the
#' \code{Splinets}--class. In the case when all the properties are satisfied the logical \code{TRUE} is returned. Otherwise,
#' \code{FALSE} is returned together with suggested corrections.
#' @param object \code{Splinets} object, the object to be diagnosed; For this object to be corrected properly each support interval has to have at least \code{2*smorder+4} knots.
#' @return A list made of: a logical value \code{is}, a \code{Splinets} object \code{robject}, a numeric value \code{Er}, and a string \code{Report}.
#' \itemize{
#' \item The logical value \code{is} indicates if all the condtions for the elements of \code{Splinets} object to be a collection of valid
#' splines are satisfied, additional diagnostic messages are printed out.
#' \item The object \code{robject} is a modified input object that has all SLOT fields modified so the conditions/restrictions
#' to be a proper spline are satisfied.
#' \item The numeric value \code{Er} is giving the total squared error of deviation of the input matrix of derivative from the conditions required for a spline. 
#' \item The string \code{Report} contains a detailed report on the diagnostic, use \code{cat(is.splinet$Report)} to print out the report. 
#' }
#' @export
#' 
#' @inheritSection Splinets-class References
#' 
#' @seealso 
#' \code{\link{Splinets-class}} for the definition of the \code{Splinets}-class;
#' \code{\link{construct}} for constructing such an object from the class;
#' \code{\link{gather}} and \code{\link{subsample}}  for combining and subsampling \code{Splinets}-objects, respectively;
#' \code{\link{plot,Splinets-method}} for plotting \code{Splinets} objects;
#'
#' @example R/Examples/ExIsSplinets.R
#'
#'
setGeneric("is.splinets",     #This creates a generic call that will be later specified by a method that depends (dispatches)
           #on the arguments of the object, see 'setMethod'
           function(object){
             standardGeneric("is.splinets")
           }
)
#' @title Diagnostics of splines
#' @description This short information is added to satisfy an R-package building requirement, see \code{\link{is.splinets}} 
#' for the actual information.  
#' @param object \code{Splinets} object, the object to be diagnosed; 
#' @rdname is.splinets-methods
#' @aliases is.splinets,Splinets-method
#' @return A list made of: a logical value \code{is}, a \code{Splinets} object \code{robject}, a numeric value \code{Er}, and a string \code{Report}..
#' \itemize{
#' \item The logical value \code{is} indicates if all the condtions for the elements of \code{Splinets} object to be a collection of valid
#' splines are satisfied, additional diagnostic messages are printed out.
#' \item The object \code{robject} is a modified input object that has all SLOT fields modified so the conditions/restrictions
#' to be a proper spline are satisfied.
#' \item The numeric value \code{Er} is giving the total squared error of deviation of the input matrix of derivative from the conditions required for a spline. 
#' \item The string \code{Report} contains a detailed report on the diagnostic, use \code{cat(is.splinet$Report)} to print out the report. 
#' }
#' 
setMethod(
  "is.splinets","Splinets", #Here the class specification of the generic "is.splinets" is made
  function(object) {
    Report="\n\nDIAGNOSTIC CHECK of a SPLINETS object"
    
    ######### KNOTS
    
    Report=paste(Report,"\n\nTHE KNOTS:  \n")
    is=TRUE  #The flag to be returned in the list
    robject=object #The modified input object that will be returned from the function

    #First: checking if there is enough knots
    n <- length(object@knots)-2
    k=object@smorder
    if (n < k) {
      Report=paste(Report,"SLOT 'knots' in a 'splinets' object should be a vector of at least length ", k+2, "\n")
      stop("Reconsider a vector of increasing values for SLOT 'knots'.\n See 'Report' field in the output list, for more information.")
    }

    #Second: checking if the order of knots is strictly increasing
    m=min(diff(object@knots))
    if(m<=0){
      is=FALSE
      robject@knots=sort(unique(object@knots))
      Report=paste(Report,"Knots are not in the strictly increasing order, which is required.\n
Ordered  knots with removed ties are given in the output `Splinets' object.\n")
    }

    #The condition for equidistant knots are checked
    eqd=diff(object@knots)
    errkn=max(eqd)-min(eqd)
    if(errkn>(0.01*object@epsilon)){
      if(robject@equid==TRUE){
        Report=paste(Report,"Absolute deviation from the equidistant state is ",errkn*100, "%\n")
        Report=paste(Report,"SLOT 'equd'  in the output splinets object is set to FALSE.\n")
      }
      robject@equid=FALSE
    }else{
      if(robject@equid==FALSE){
        Report=paste(Report,"Absolute deviation from the equidistant state is ",errkn*100, "%.\n")
        Report=paste(Report,"SLOT 'equd'  in the output splinets object is set to TRUE.\n")
      }
      robject@equid=TRUE
    }

    #Evaluation of the matrix used in the Taylor expansion at the knots
    if(min(dim(as.matrix(object@taylor)))==1){
      if(k!=0){             #one dimensional vector has to be treated as a column matrix
        object@taylor=matrix(object@taylor,nrow=1)
      }else{ #one dimensional vector has to be treated as a row matrix
        object@taylor=matrix(object@taylor,ncol=1)
      }

    }#To assure that illogical treatment of vectors and matices in R is fixed.


    if(dim(object@taylor)[2]!=(k+1)){ #If the matrix has been evaluated before, then it has 'k+1' columns
      Report=paste(Report,"The Taylor expansion coefficient matrix does not have the proper number of columns.\n")
      is=FALSE
      if(robject@equid==TRUE){
        robject@taylor<-taylor_coeff(object@knots[1:2],k) #all the coefficients for all other knots are the same as for
        #the first two
      }else{
        robject@taylor<-taylor_coeff(object@knots,k)
      }
      Report=paste(Report,"The Taylor coefficient has been evaluated and assigned to the output 'Splinets-object'.\n")
    }else{
      if(robject@equid==FALSE){
        if(dim(object@taylor)[1]!=(n+1)){#If the matrix has been evaluated before, then it has 'n+1' columns in the non-equidistant case
          Report=paste(Report,"The proper Taylor expansion coefficient matrix does not have the proper number of rows.\n")
          is=FALSE
          robject@taylor<-taylor_coeff(object@knots,k)
          Report=paste(Report,"The Taylor coefficient has been evaluated and assigned to the output.\n")
          }
      }else{
        if(dim(object@taylor)[1]!=1){#If the matrix has been evaluated before, then it has 'n+1' columns in the non-equidistant case
          Report=paste(Report,"The proper Taylor expansion coefficient matrix does not have the proper number of rows.\n")
          is=FALSE
          robject@taylor<-taylor_coeff(object@knots[1:2],k) #all the coefficients for all other knots are the same as for
          #the first two
          Report=paste(Report,"The Taylor coefficient has been evaluated and assigned to the output.\n")
        }
      }
    }
    
    ################ SUPPORT SETS

    Report=paste(Report,"\n\nTHE SUPPORT SETS:  \n\n")
    
    #Two cases: full support for all splines and not full support
    
    sz=length(object@supp) 
 
    if(sz==0){ #CASE 1:the full support
      
      Report=paste(Report,"The support sets for the splines are equal to the entire range of knots.\n")
      robject@supp=object@supp # the empty list in this case
      
      #END CASE 1
      }else{  #CASE 2: different not full support sets for some splines
        
      for(i in 1:sz){ #running through all splines as 'sz' is equal to the number of splines in the splinet
        ss=object@supp[[i]] #the support set for the current spline
        if(min(dim(t(ss)))==1){ss=matrix(ss,ncol=2)} #the case of a single interval to be treated as a row matrix
        Nsup=dim(ss)[1] #The number of support intervals in the support set
        le=0 #to mark ending of the preceeding interval in the loop
        for(j in 1:Nsup) {#Running through the intervals of the support
          if((ss[j,1]-le)<1){
            Report=paste(Report,"SLOT 'supp' in the input 'splinets' object is not valid.\n The rows in the matrices should form separated disjoint intervals but they overlap.\n")
            stop("Review the values in SLOT 'supp'.\n See 'Report' field in the output list, for more information.")
          }
          le=ss[j,2]
          if((ss[j,2]-ss[j,1])<(k+1)){
            Report=paste(Report,"SLOT 'supp' in the input 'splinets' object is not valid. \n The rows in the matrices should correspond to", k+2, "knots in each support interval.\n")
            stop("Not enough knots per interval in the support, review the values in SLOT 'supp'.\n See 'Report' field in the output list, for more information.")
          }
        }
        robject@supp[[i]]=ss 
        
      #END CASE 2
      }
    }
    ##################The most important fields of the splinets class###############

    Report=paste(Report,"\n\nTHE DERIVATIVES AT THE KNOTS:  \n\n")
        
    #The diagnostic has two parts: Part I: Checking the dimentions of the matrices of derivatives at the knots
                                #  Part II: Checking the boundary conditions for the matrices
                                #  Part III: Checking the main conditions for the matrices
        
    #PART I - the dimensions of the matrices
    #There are two cases full support and not full
    size=length(object@der)
    
    if(sz==0){ #CASE 1  full support for all splines
      for(j in 1:size){#Running through all the splines in the object
        if(prod(dim(object@der[[j]])==c(n+2,k+1))==0){
          is=FALSE
          Report=paste(Report,"The full support range case.\n")
          Report=paste(Report,"SLOT 'der' do not have properly set the dimension for spline", j,"in the input 'Splinets' object.\n")
          Report=paste(Report,"In the output object, it is set temporarily to the matrix of 'ones' with the proper dimension.\n\n")
          robject@der[[j]]=matrix(1,nrow=n+2,ncol=k+1)
        }
      }
      #END: CASE 1
    }else{#CASE 2 arbitrary supports
      if(sz!=size){stop("The sizes of SLOTS 'supp' and 'der' should match, the object needs correction.\n See 'Report' field in the output list, for more information.")}
      for(j in 1:size){#Running through all the splines in the object
          #In the condition on the RHS below the sum of all knots in the support set is evaluated and should match the number
          #of rows in 'der'.
          m=sum(object@supp[[j]][,2]-object@supp[[j]][,1]+1) #The total number of knots in the support
        if(prod(dim(object@der[[j]])==c(m,k+1))==0){
          is=FALSE
          Report=paste(Report,"The partial support range case.\n")
          Report=paste(Report,"SLOT 'der' do not have properly set the dimension for spline", j,"in the input 'Splinets' object.\n")
          Report=paste(Report,"In the output object, it is set temporarily to the matrix of 'ones' with the proper dimension.\n\n")
          robject@der[[j]]=matrix(1,nrow=m,ncol=k+1)
        }
      }
     }
    #End of PART I - the dimensions of the matrices    
  
  Er=0 #the initial value for the deviation from the condtion   
  
  if(k>0){ #There are no boundary or other condtions in the zero order case. 
    
   #PART II - the boundary conditions for the matrices
        #There are two cases to consider: the full support and non-full support
  
      
        if(sz==0){ #CASE 1 full support for all splines
            for(j in 1:size){#Running through all the splines in the object
            CurrDer=robject@der[[j]]
            if(prod(abs(CurrDer[1,(1:k)])<rep(object@epsilon,k))==0 || prod(abs(CurrDer[dim(CurrDer)[1],(1:k)])<rep(object@epsilon,k))==0){
              is=FALSE
              Report=paste(Report,"The boundary zero conditions are not satisfied for spline", j , "in the input 'Splinets' object.\n")
              robject@der[[j]][1,(1:k)]=rep(0,k)
              robject@der[[j]][dim(CurrDer)[1],(1:k)]=rep(0,k)
              Report=paste(Report,"Correction of the first and last rows of the derivative matrices are made in the output 'Splinets' object.\n")
            }
          } #The end of verification of the boundary conditions
          #END: CASE 1
        }else{#CASE 2 non-full support for all splines
          for(j in 1:size){#Running through all the splines in the object
            CurrDer=robject@der[[j]]
            CurrSupp=object@supp[[j]]
            B=1 #Counting cumulatively knots in the support
            Nsupp=dim(CurrSupp)[1]
            for(l in 1:Nsupp){#running through the support intervals
              E=CurrSupp[l,2]-CurrSupp[l,1]+B
            if(prod(abs(CurrDer[B,(1:k)])<rep(object@epsilon,k))==0 || prod(abs(CurrDer[E,(1:k)])<rep(object@epsilon,k))==0){
              is=FALSE
              Report=paste(Report,"The boundary zero conditions are not satisfied for spline", j , "in the input 'Splinets' object.\n")
              robject@der[[j]][B,(1:k)]=rep(0,k)
              robject@der[[j]][E,(1:k)]=rep(0,k)
              Report=paste(Report,"Correction of the first and last rows of the derivative matrices over the support component", l,"of spline", j,"in the output 'Splinets' object.\n")
            }
              B=E+1
           }
          } #The end of verification of the boundary derivative matrices
          #END: CASE 2
        }

        #End of PART II - the boundary conditions for the matrices
    
    #PART III - the main conditions for the matrices
    
    #There are four cases to consider controlled by two factors at two levels: the equidistant and full support
    
    if(sz==0 && object@equid==TRUE){ #CASE 1 equidistant and full support for all splines
      l=floor(n/2)
      #Evaluating the main condition for the proper spline object:
      AL=robject@taylor 
      AR=AL%*%diag((-1)^(0:k))
      
      for(j in 1:size){
      
      SL=as.matrix(robject@der[[j]][1:(l+2),]) #These matrices are explained in the paper. They go one row beyond
      SR=as.matrix(robject@der[[j]][(n+2)-(0:(l+1)),]) #the midpoint in the even case, in the odd case they share the center
                                                       #row of the 'der' matrix as their last rows.
      
      #Part 1a: The higest derivative at the center of the derivative matrices ###OBS!!!
      
      if(n %% 2 == 0){ #The even number of knots
        if(SL[l+2,k+1]!=SR[l+2,k+1]){
          Report=paste(Report,"\nSpline", j,"'s highest derivative is not symmetrically defined at the center (the values at the two central knots should be equal).\n")
          Report=paste(Report,"Spline", j,"'s highest derivative values at the two central knots have been made equal by averaging the two central values in SLOT 'der'.\n")
          is=FALSE
          robject@der[[j]][(l+2),k+1]=(SL[l+2,k+1]+SR[l+2,k+1])/2  #Correcting values at the center
          robject@der[[j]][(l+1),k+1]=robject@der[[j]][(l+2),k+1]
          SL[l+2,k+1]=robject@der[[j]][(l+2),k+1]
          SR[l+2,k+1]=SL[l+2,k+1]
        }
      }else{ #The odd number of knots
        if((SL[l+2,k+1])^2+(SR[l+2,k+1])^2!=0){
          Report=paste(Report,"\nSpline", j,"'s highest derivative at the central knot is not equal to zero.\n")
          Report=paste(Report,"Spline", j,"'s highest derivative value at the central knot has been made equal to zero.\n")
          is=FALSE
          SL[l+2,k+1]=0
          SR[l+2,k+1]=0
          robject@der[[j]][(l+2),(k+1)]=0
        }
      }
        
        Er1=derver(SL,AL)+derver(SR,AR)
        Er=Er+Er1 #Adding the error for the j-th spline to the total error over all splines
        #The main verification of the conditions and the total squared error of the deviation.
        
        if(sqrt(Er1/(n*k))>object@epsilon){
          is=FALSE
          Report=paste(Report,"\nThe matrix of derivatives at the knots for spline", j, "does not satisfy the conditions \n
          required for a spline (up to the accuracy SLOT 'epsilon').\n
          One of the reasons can be that SLOT 'taylor' is not correctly given.\n")
          Report=paste(Report,"The computed standard error per matrix entry is", sqrt(Er1/(n*k)), ".\n\n")
          #Recomputing the matrix of the derivative using RRM method 
            robject@der[[j]]=correct(robject@der[[j]],robject@taylor)  
            Report=paste(Report,"\nThe output object Spline", j," has the derivative matrix corrected by the RRM method\n given that SLOT 'taylor' is properly given.")
        }
      }
      #END: CASE 1
    }else{
      if(object@equid==TRUE){ #CASE 2 eqiudistant and arbitrary supports
       
        #Evaluating the main condition for the proper spline object:
        AL=robject@taylor 
        AR=AL%*%diag((-1)^(0:k))
        
        for(j in 1:size){
          CurrDer=robject@der[[j]]
          CurrSupp=object@supp[[j]]
          B=1 #Counting cumulatively knots in the support
          Nsupp=dim(CurrSupp)[1]
          for(l in 1:Nsupp){#running through the support intervals
          m=CurrSupp[l,2]-CurrSupp[l,1]-1 #the number of internal knots in the l-th interval of the support
          L=floor(m/2)
          E=CurrSupp[l,2]-CurrSupp[l,1]+B  #the index for the entries in 'der' matrix for the last knot in the current support interval
          SL=as.matrix(robject@der[[j]][B:(B+L+1),]) #These matrices are explained in the paper. They go one row beyond
          SR=as.matrix(robject@der[[j]][E-(0:(L+1)),]) #the midpoint in the even case, in the odd case they share the center
          #row of the 'der' matrix as their last rows.
          
          #Part 1a: The higest derivative at the center of the derivative matrices ###OBS!!!
          
          if(m %% 2 == 0){ #The even number of knots
            if(SL[L+2,k+1]!=SR[L+2,k+1]){
              Report=paste(Report,"\nSpline", j,", support", l,"'s highest derivative is not symmetrically defined at the center (the values at the two central knots should be equal).\n")
              Report=paste(Report,"Spline", j," highest, support", l,"'s derivative values at the two central knots have been made equal by averaging the two central values in SLOT 'der'.\n")
              is=FALSE
              robject@der[[j]][(B+L+1),k+1]=(SL[L+2,k+1]+SR[L+2,k+1])/2  #Correcting values at the center
              robject@der[[j]][(B+L),k+1]=robject@der[[j]][(B+L+1),k+1]
              SL[L+2,k+1]=robject@der[[j]][(B+L+1),k+1]
              SR[L+2,k+1]=SL[L+2,k+1]
            }
          }else{ #The odd number of knots
            if((SL[L+2,k+1])^2+(SR[L+2,k+1])^2!=0){
              Report=paste(Report,"\nSpline", j,"support", l,"'s highest derivative at the central knot is not equal to zero.\n")
              Report=paste(Report,"Spline", j,"support", l,"'s highest derivative value at the central knot has been made equal to zero.\n")
              is=FALSE
              SL[L+2,k+1]=0
              SR[L+2,k+1]=0
              robject@der[[j]][(B+L+1),(k+1)]=0
            }
          }
          
          Er1=derver(SL,AL)+derver(SR,AR) #The main verification of the conditions and the standard devation from them per the matrix entry.
          Er=Er+Er1 #Adding the error for the j-th spline to the total error over all splines
          if(sqrt(Er1/(m*k))>object@epsilon){
            is=FALSE
            Report=paste(Report,"\nThe matrix of derivatives at the knots for spline", j, ", support", l," does not satisfy the conditions that are required for a spline (up to the accuracy SLOT 'epsilon').\n")
            Report=paste(Report,"The computed standard error per matrix entry is", sqrt(Er1/(m*k)), ".\n\n")
            #Recomputing the matrix of the derivative using RRM method 
            robject@der[[j]][B:E,]=correct(robject@der[[j]][B:E,],robject@taylor)  
            Report=paste(Report,"\nThe output object Spline", j," support", l,"has the derivative matrix corrected by the RRM method.")
          }
          B=E+1
          }#the end of the loop over the support intervals
        }#the end of the loop over the splines
        #END: CASE 2
      }else{
        if(sz==0){#CASE 3 non-equidistant and full support
          l=floor(n/2)
          for(j in 1:size){
            SL=as.matrix(robject@der[[j]][1:(l+2),]) #These matrices are explained in the paper. They go one row beyond
            SR=as.matrix(robject@der[[j]][(n+2)-(0:(l+1)),]) #the midpoint in the even case, in the odd case they share the center
            #row of the 'der' matrix as their last rows.
            
            #Part 1a: The higest derivative at the center of the derivative matrices ###OBS!!!
            
            if(n %% 2 == 0){ #The even number of knots
              if(SL[l+2,k+1]!=SR[l+2,k+1]){
                Report=paste(Report,"\nSpline", j,"'s highest derivative is not symmetrically defined at the center (the values at the two central knots should be equal).\n")
                Report=paste(Report,"The spline", j,"'ths highest derivative at the two central knots has been made equal by averaging SLOT 'der'.\n")
                is=FALSE
                robject@der[[j]][(l+2),k+1]=(SL[l+2,k+1]+SR[l+2,k+1])/2  #Correcting values at the center
                robject@der[[j]][(l+1),k+1]=robject@der[[j]][(l+2),k+1]
                SL[l+2,k+1]=robject@der[[j]][(l+2),k+1]
                SR[l+2,k+1]=SL[l+2,k+1]
              }
            }else{ #The odd number of knots
              if((SL[l+2,k+1])^2+(SR[l+2,k+1])^2!=0){
                Report=paste(Report,"\n The spline", j,"'ths highest derivative at the central knot is zero.\n")
                Report=paste(Report,"Now it is set to zero.\n")
                is=FALSE
                SL[l+2,k+1]=0
                SR[l+2,k+1]=0
                robject@der[[j]][(l+2),(k+1)]=0
              }
            }
            #Evaluating the main condition for the proper spline object:
            AL=as.matrix(robject@taylor[(1:(l+1)),]) #The top part of the split of 'taylor' matrix
            AR=as.matrix(robject@taylor[((n+1):(n-l+1)),]%*%diag((-1)^(0:k))) #The bottom part of the split in the reverse order
            #of 'taylor' matrix
            #Matrices 'AL' and 'AR' will be used for verification of the consistency of
            #the conditions for the derivatives
            
            Er1=derver(SL,AL)+derver(SR,AR) #The main verification of the conditions and the standard devation from them per the matrix entry.
            Er=Er+Er1 #Adding the error for the j-th spline to the total error over all splines
            if(sqrt(Er1/(n*k))>object@epsilon){
              is=FALSE
              Report=paste(Report,"\nThe derivative matrix for spline", j, "does not satisfy the smoothness conditions (up to the accuracy SLOT 'epsilon').\n")
              Report=paste(Report,"The standard error per matrix entry is", sqrt(Er1/(n*k)), ".\n\n")
              #Recomputing the matrix of the derivative using RRM method 
              robject@der[[j]]=correct(robject@der[[j]],robject@taylor)  
              Report=paste(Report,"\nThe output object has the derivative matrix corrected by the RRM method.\n")
            }
          }  
          #END: CASE 3
        }else{#CASE 4 non-equidistant and arbitrary supports
          #Evaluating the main condition for the proper spline object:
          for(j in 1:size){
            CurrDer=robject@der[[j]]
            CurrSupp=object@supp[[j]]
            B=1 #Counting cumulatively knots in the support
            Nsupp=dim(CurrSupp)[1]
            for(l in 1:Nsupp){#running through the support intervals
              m=CurrSupp[l,2]-CurrSupp[l,1]-1 #the number of internal knots in the l-th interval of the support
              L=floor(m/2)
              E=CurrSupp[l,2]-CurrSupp[l,1]+B  #the index for the entries in 'der' matrix for the last knot in the current support interval
              SL=as.matrix(robject@der[[j]][B:(B+L+1),]) #These matrices are explained in the paper. They go one row beyond
              SR=as.matrix(robject@der[[j]][E-(0:(L+1)),]) #the midpoint in the even case, in the odd case they share the center
              #row of the 'der' matrix as their last rows.
              
              #Part 1a: The higest derivative at the center of the derivative matrices ###OBS!!!
              
              if(m %% 2 == 0){ #The even number of knots
                if(SL[L+2,k+1]!=SR[L+2,k+1]){
                  Report=paste(Report,"\nSpline", j,", support", l," - highest derivative is not symmetric at the center (equal values at the two central knots).\n")
                  Report=paste(Report,"The two values have been made equal by averaging.\n")
                  is=FALSE
                  robject@der[[j]][(B+L+1),k+1]=(SL[L+2,k+1]+SR[L+2,k+1])/2  #Correcting values at the center
                  robject@der[[j]][(B+L),k+1]=robject@der[[j]][(B+L+1),k+1]
                  SL[L+2,k+1]=robject@der[[j]][(B+L+1),k+1]
                  SR[L+2,k+1]=SL[L+2,k+1]
                }
              }else{ #The odd number of knots
                if((SL[L+2,k+1])^2+(SR[L+2,k+1])^2!=0){
                  Report=paste(Report,"\nSpline", j,"support", l,"'s highest derivative at the central knot is not zero.\n")
                  Report=paste(Report,"Now it is set to zero.\n")
                  is=FALSE
                  SL[L+2,k+1]=0
                  SR[L+2,k+1]=0
                  robject@der[[j]][(B+L+1),(k+1)]=0
                }
              }
              #Evaluating the main condition for the proper spline object:
              AL=as.matrix(robject@taylor[CurrSupp[l,1]:(CurrSupp[l,1]+L),]) #The top part of the split of 'taylor' matrix
              AR=as.matrix(robject@taylor[(CurrSupp[l,2]-1):(CurrSupp[l,2]-L-1),]%*%diag((-1)^(0:k))) #The bottom part of the split in the reverse order
              #of 'taylor' matrix
              #Matrices 'AL' and 'AR' will be used for verification of the consistency of
              #the conditions for the derivatives
              
              Er1=derver(SL,AL)+derver(SR,AR) #The main verification of the conditions and the standard devation from them per the matrix entry.
              Er=Er+Er1 #Adding the error for the j-th spline to the total error over all splines
              if(sqrt(Er1/(m*k))>object@epsilon){
                is=FALSE
                Report=paste(Report,"\nThe matrix of derivatives at the knots for spline", j, ", support", l," does not satisfy the splie conditions (up to the accuracy set in SLOT 'epsilon').\n")
                cat("The computed standard error per matrix entry is", sqrt(Er1/(m*k)), ".\n\n")
                #Recomputing the matrix of the derivative using RRM method 
                robject@der[[j]][B:E,]=correct(robject@der[[j]][B:E,],robject@taylor[CurrSupp[l,1]:(CurrSupp[l,2]-1),])  
                Report=paste(Report,"\nThe output object Spline", j," support", l,"has the derivative matrix corrected by the RRM method.")
              }
              B=E+1
            }#the end of the loop over the support intervals
          }#the end of the loop over the splines
          #END: CASE 4
        }
      }
    }
  } #The end of the non-zero smoothness order.
    #End of PART III - the main condtions for the matrices

    issp=list(is,robject,Er,Report)
    names(issp)=c("is","robject","Er","Report")
    if(is==TRUE){
      message("\n The input object is a 'Splinets' object.")
    }else{
        message("\n The input object is not a 'Splinets' object.\n See 'Report' field in the output list, for more information.")
      }
    return(invisible(issp))     

 } #The end of the method function
) #The end of the method
