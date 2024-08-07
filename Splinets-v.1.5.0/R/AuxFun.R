#--- Functions are grouped with respect to the part of the package they are used ---#
#-------------------- Within the group the order is alphabetic ---------------------#

#------------------------------#
#--- Sec 1 AuxFun for Class ---#
#------------------------------#

# For a given matrix returns the corrected matrix that corresponds to a spline using one of 
# the three methods described in the paper.
# S -- (n+2) x (k+1) matrix. The matrix in the symmetric around the center to be corrected according to one of
#     the three methods described in the paper.
#      n>=2*k.
# Vd -- (n+1) x (k+1) matrix if non-equidistant knots, or 1 x (k+1) if the equidistant knots (vectors are fine)
#      columnwise vectors  of the Taylor expansion coefficients at the knots,
# method -- string. One of the following: 'CRLC', 'CRFC', 'RRM' (Default), indicating the method for correcting the input matrix S.
#
# return a matrix (m+2) x (k+1) SS filled with entries of the derivatives at the knots, in the symmetric format
# and with the zero boundary condtions, that satisfies the restrictions.
correct=function(S,Vd,method='RRM'){
  k=dim(S)[2]-1
  n=dim(S)[1]-2
  l=floor(n/2)
  SS=S

  if(n<(2*k+2)){
  stop("\nThe number of the support knots is too small for the considered matrix correction methods.\n Consider increasing their number to be at least ",2*k+4,'.\n')
    }
  if(min(dim(as.matrix(Vd)))==1){Vd=matrix(Vd,nrow=1)} #To treat vector inputs of Vd (fixing R treatment of vectors and matrices)
  if(dim(Vd)[2]!=(k+1)){stop("Incompatible dimensions of the input matrices.")}
  if(dim(Vd)[1]>1){
    if(dim(Vd)[1]!=(n+1)){stop("Incompatible dimensions of the input matrices.")}
  }

  #Checking boundary conditions
  if(k>0){ #There are no boundary condtions in the zero order case.
    if(prod(S[1,(1:k)]==rep(0,k))==0 || prod(S[dim(S)[1],(1:k)]==rep(0,k))==0){
      cat("The zero boundary conditions are not satisfied.\n")
      SS[1,(1:k)]=rep(0,k)
      SS[dim(SS)[1],(1:k)]=rep(0,k)
      cat("The correction of the first and last rows of the derivative matrix has been made.\n\n")
    }
  } #The end of verification of the boundary derivative matrices


  #Verification of the conditions at the center of the derivative matrices
  if(n %% 2 == 0){
    if(S[l+2,k+1]!=S[l+1,k+1]){
      cat("\nThe highest order derivative is not symmetrically defined at the center.
          The values at the two central knots should be equal.\n")
      cat("The highest order derivative values at the two central knots
          have been made equal by averaging the two central values.\n")
      SS[(l+2),k+1]=(S[l+2,k+1]+S[l+1,k+1])/2  #Correcting values at the center
      SS[l+1,k+1]=SS[(l+2),k+1]
    }
  }else{
    if((S[l+2,k+1])^2!=0){
      cat("\nThe highest order derivative at the central knot is not equal to zero.\n")
      cat("It has been made equal to zero now.\n")
      SS[l+2,(k+1)]=0
    }
  }


  #Next we correct the two central rows so they correspond to a valid spline over the central interval
  #(the case of n even only) using 'frlr' routine with m=0 (two rows) and averaging the results

  if(n %% 2 == 0){
    if(dim(Vd)[1]==1){
      Vc=Vd
    }else{
      Vc=Vd[l+1,]
    }
    UU=frlr(SS[l+1,],SS[l+2,],Vc)
    VV=frlr(SS[l+2,],SS[l+1,],Vc%*%diag((-1)^(0:k)))
    UU[1,]=(VV[2,]+UU[1,])/2
    UU[2,]=(VV[1,]+UU[2,])/2
    SS[(l+1):(l+2),]=UU
  }
  #At this point matrix SS has the proper boundary conditions and has well defined rows at the center.
  #The rest of the entries needs to be fixed next through the split the matrix to halves.

  #Dividing the matrix into the left and right halves
  SL=as.matrix(SS[1:(l+2),])
  SR=as.matrix(SS[(n+2):(n-l+1),])

  #Definition of the left and the right Taylor matrices
  if(dim(Vd)[1]!=1){
    AL=as.matrix(Vd[1:(l+1),]) #The top part of the split of 'taylor' matrix

    AR=as.matrix(Vd[(n+1):(n+1-l),]%*%diag((-1)^(0:k))) #The bottom part of the split

    #of 'taylor' matrix
  }else{
    AL=Vd
    AR=AL%*%diag((-1)^(0:k))
  }

  #Since the correction of the matrices will go from the center to the endpoints, we reverse the order of columns.

  SL=SL[(l+2):1,]
  SR=SR[(l+2):1,]

  if(dim(Vd)[1]!=1){ #Equally spaced vs regular
    AL=AL[(l+1):1,]%*%diag((-1)^(0:k))
    AR=AR[(l+1):1,]%*%diag((-1)^(0:k))
  }else{
    AL=AL%*%diag((-1)^(0:k))
    AR=AR%*%diag((-1)^(0:k))
  }

  AAL=AL #Reversed commplete Taylor coefficient matrices for each half
  AAR=AR



  if(method=='CRLC'){  #The second method, the center row, the last column
    if((dim(AL)[1])!=1){
      AL=AL[1:(l-k),]
    }
    if((dim(AR)[1])!=1){
      AR=AR[1:(l-k),]
    }

    FR=SL[1,]
    LC=SL[2:(l-k+2),k+1]
    FR[k+1]=LC[1]
    SL[1:(l-k+1),]=frlc(FR,LC,AL)

    FR=SR[1,]
    LC=SR[2:(l-k+2),k+1] #The necessary shift by 1 of the highest order derivative values due to the reversed order
    FR[k+1]=LC[1]         #the change of the direction in the matrices
    SR[1:(l-k+1),]=frlc(FR,LC,AR)

    #Now it only left to evaluate the derivatives at the final k+1 knots at both ends, using 'frlr' method.
    if(dim(AL)[1]!=1){
      AL=AAL[(l-k+1):(l+1),]
    }
    if(dim(AR)[1]!=1){
      AR=AAR[(l-k+1):(l+1),]
    }
    
    cat('\nCorrection of the LHS part of the matrix')
    SL[(l-k+1):(l+2),]=frlr(SL[(l-k+1),],SL[(l+2),],AL,neqknots=k)
    cat('\nCorrection of the RHS part of the matrix')
    SR[(l-k+1):(l+2),]=frlr(SR[(l-k+1),],SR[(l+2),],AR,neqknots=k)

    #Puting the output from the computed pieces
    if(n %% 2 == 0 ){
      SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
      SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]

      SS[(l+1):(2*l+2),1:k]=SR[1:(l+2),1:k]
      SS[(l+2):(2*l+2),(k+1)]=SR[1:(l+1),(k+1)]
    }else{
      SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
      SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]
      SS[l+2,k+1]=0
      SS[(l+2):(2*l+3),1:k]=SR[1:(l+2),1:k]
      SS[(l+3):(2*l+3),(k+1)]=SR[1:(l+1),(k+1)]
    }

  }else{                #The third method center row, first column
    if(method=='CRFC'){
      if((dim(AL)[1])!=1){
        AL=AL[1:(l-k+1),,drop = FALSE]
      }
      if((dim(AR)[1])!=1){
        AR=AR[1:(l-k+1),,drop = FALSE]
      }

      FR=SL[1,1:k]
      FC=SL[1:(l-k+2),1]

      SL[1:(l-k+2),]=frfc(FR,FC,AL)

      FR=SR[1,1:k]
      FC=SR[1:(l-k+2),1]

      SR[1:(l-k+2),]=frfc(FR,FC,AR)

      #Now it only left to evaluate the derivatives at the final k+1 knots at both ends, using 'frlr' method.
      if(dim(AL)[1]!=1){
        AL=AAL[(l-k+1):(l+1),,drop = FALSE]
      }
      if(dim(AR)[1]!=1){
        AR=AAR[(l-k+1):(l+1),,drop = FALSE]
      }
      cat('\nCorrection of the LHS part of the matrix')
      SL[(l-k+1):(l+2),]=frlr(SL[(l-k+1),],SL[(l+2),],AL,neqknots=k)
      cat('\nCorrection of the RHS part of the matrix')
      SR[(l-k+1):(l+2),]=frlr(SR[(l-k+1),],SR[(l+2),],AR,neqknots=k)

      #Puting the output from the computed pieces
      if(n %% 2 == 0 ){
        SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
        SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]

        SS[(l+1):(2*l+2),1:k]=SR[1:(l+2),1:k]
        SS[(l+2):(2*l+2),(k+1)]=SR[1:(l+1),(k+1)]
      }else{
        SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
        SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]
        SS[l+2,k+1]=0
        SS[(l+2):(2*l+3),1:k]=SR[1:(l+2),1:k]
        SS[(l+3):(2*l+3),(k+1)]=SR[1:(l+1),(k+1)]
      }
    }else{ # The first method, or the regular row method (RRM) - the default

      j=floor((l-k)/(k+1)) #This determines the number blocks for `frlr` with m=k


      SL[1:(l+1),k+1]=SL[2:(l+2),k+1] #the shift of the vector of the highest order derivative
      #by one due to reversion of the direction of knots
      SR[1:(l+1),k+1]=SR[2:(l+2),k+1]

      if(j>0){#only if there is at least one full block with k+2 knots between the center and the final one

        for(r in 1:j){   #The loop of full blocks of the sizes k+2
          FR=SL[(1+(r-1)*(k+1)),]
          LR=SL[(1+r*(k+1)),]
          if((dim(AAL)[1])!=1){
            AL=AAL[(1+(r-1)*(k+1)):(r*(k+1)),]
          }
          cat('\nCorrection of the LHS part of the matrix')
          SL[(1+(r-1)*(k+1)):(1+r*(k+1)),]=frlr(FR,LR,AL,neqknots=k)

          FR=SR[(1+(r-1)*(k+1)),]
          LR=SR[(1+r*(k+1)),]
          if((dim(AAR)[1])!=1){
            AR=AAR[(1+(r-1)*(k+1)):(r*(k+1)),]
          }
          cat('\nCorrection of the RHS part of the matrix')
          SR[(1+(r-1)*(k+1)):(1+r*(k+1)),]=frlr(FR,LR,AR,neqknots=k)

        }

      }

      if((j*(k+1)+1) != (l-k+1)){ #the case of one incomplete block of knots to be handled by 'frlr' with m<k

        FR=SL[j*(k+1)+1,]
        LR=SL[l-k+1,]
        if((dim(AAL)[1])!=1){
          AL=AAL[(j*(k+1)+1):(l-k),]
        }
        cat('\nCorrection of the LHS part of the matrix')
        SL[(j*(k+1)+1):(l-k+1),]=frlr(FR,LR,AL,neqknots=(l-(k+1)*(j+1)))

        FR=SR[j*(k+1)+1,]
        LR=SR[l-k+1,]
        if((dim(AAR)[1])!=1){
          AR=AAR[(j*(k+1)+1):(l-k),]
        }
        cat('\nCorrection of the RHS part of the matrix')
        SR[(j*(k+1)+1):(l-k+1),]=frlr(FR,LR,AR,neqknots=(l-(k+1)*(j+1)))

      }

      #Now it only left to evaluate the derivatives at the final k+1 knots at both ends, using 'frlr' method.
      if(dim(AAL)[1]!=1){
        AL=AAL[(l-k+1):(l+1),]
      }
      if(dim(AAR)[1]!=1){
        AR=AAR[(l-k+1):(l+1),]
      }
      cat('\nCorrection of the LHS part of the matrix')
      SL[(l-k+1):(l+2),]=frlr(SL[(l-k+1),],SL[(l+2),],AL,neqknots=k)
      cat('\nCorrection of the RHS part of the matrix')
      SR[(l-k+1):(l+2),]=frlr(SR[(l-k+1),],SR[(l+2),],AR,neqknots=k)

      #Puting the output from the computed pieces
      if(n %% 2 == 0 ){
        SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
        SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]

        SS[(l+1):(2*l+2),1:k]=SR[1:(l+2),1:k]
        SS[(l+2):(2*l+2),(k+1)]=SR[1:(l+1),(k+1)]
      }else{
        SS[1:(l+2),1:k]=SL[(l+2):1,1:k]
        SS[1:(l+1),(k+1)]=SL[(l+1):1,(k+1)]
        SS[l+2,k+1]=0
        SS[(l+2):(2*l+3),1:k]=SR[1:(l+2),1:k]
        SS[(l+3):(2*l+3),(k+1)]=SR[1:(l+1),(k+1)]
      }
    }
  }

  SS[1,1:k]=rep(0,k) #To remove computational inaccuracies at the boundaries
  SS[n+2,1:k]=rep(0,k)

  return(SS)
}

#The following function verifies the knot conditions for the matrix of the derivatives.
#S -- a (l+2) x (k+1) matrix,
#Vd -- a (l+1) x (k+1) matrix if non-equidistant knots, or 1 x (k+1) if the equidistant knots (vectors are fine)
#      columnwise vectors  of the Taylor expansion coefficients at the knots,
#Returns the mean squared error.
derver=function(S,Vd){
  Vd=as.matrix(Vd) #To assure that dim(Vd) works properly for the vector input

  l=dim(S)[1]-2
  k=dim(S)[2]-1
  y=S

  if(dim(Vd)[1]==1){
    if((k+1)!=dim(Vd)[2]){stop("The size of the input vector incompatible with the size of the input matirx.")}
    y[2:(l+2),]=S[1:(l+1),]%*%toeplitz_lower(Vd) #It does not matter if 'Vd' is a column, row or vector.
  }else{
    if(k!=dim(Vd)[2]-1){stop("Incompatible number of columns in the input matrices.")}
    if(l!=dim(Vd)[1]-1){stop("Incompatible number of rows in the input matrices.")}
      for(i in 1:(l+1)){
        y[i+1,]=S[i,]%*%toeplitz_lower(Vd[i,])
      }
  }

  y[2:(l+2),k+1]=y[2:(l+2),k+1]+diff(S[,k+1])
  Er=sum((y-S)^2)
  return(Er)
}

#First row -- first column verification as described in the paper
# FR --  k  numeric vector. The first row in the output matrix, except for the highest (k-th) derivative
# FC -- (m+2) numeric vector. The last column in the output matrix,
# Vd -- a (l+1) x (k+1) matrix if non-equidistant knots, or 1 x (k+1) if the equidistant knots (vectors are fine)
#      columnwise vectors  of the Taylor expansion coefficients at the knots,
# LRL -- an optional value of the last knot, the k-th derivative (not involved in the computations and restrictions),
#      it is set to the optional value.
# return a matrix (m+2) x (k+1) S filled with entries of the derivatives at the knots.
#WARNING! This solution should not be used for large n because it becomes very unstable (tests indicate anything above
#20 knots starts fail)
frfc=function(FR,FC,Vd,LRL=0){
  k=length(FR)
  m=length(FC)-2
  if(FR[1]!=FC[1])
  {
    stop("The first entry of the 'first row' argument is not the same as the first entry of the 'first column' argument.")
  }

  S=matrix(0,ncol=k+1,nrow=m+2)
  S[1,1:k]=FR
  S[,1]=FC
  S[m+2,k+1]=LRL #Setting the lower right corner
  Vd=as.matrix(Vd) #To assure that dim(Vd) works properly for the vector input
  if(dim(Vd)[2]==1){Vd=t(Vd)}
  if((k+1)!=dim(Vd)[2]){stop("Incompatible dimensions of the input matrices.")}

  #Evaluating all the other entries in S
  if(dim(Vd)[1]==1){ #the equally spaced knots case
    A=toeplitz_lower(Vd[1,])
    for(i in 2:(m+2)){
      S[i-1,k+1]=S[i,1]/A[k+1,1]-S[i-1,1:k]%*%(A[1:k,1]/A[k+1,1])
      PS=S[i-1,]%*%A
      S[i,2:k]=PS[2:k]
    }
  }else{
    if((m+1)!=dim(Vd)[1]){stop("Incompatible dimensions of the input matrices.")}
    for(i in 2:(m+2)){
      A=toeplitz_lower(Vd[i-1,])
      S[i-1,k+1]=S[i,1]/A[k+1,1]-S[i-1,1:k]%*%(A[1:k,1]/A[k+1,1])
      PS=S[i-1,]%*%A
      S[i,2:k]=PS[2:k]
    }
  }
  return(S)
}

#First row -- last column construction of the matrix of the derivatives at knots as described in the paper
# FR -- (k+1) numeric vector. The first row in the output matrix
# LC -- (m+2) numeric vector. The last column in the output matrix
# Vd -- a (l+1) x (k+1) matrix if non-equidistant knots, or 1 x (k+1) if the equidistant knots (vectors are fine)
#      columnwise vectors  of the Taylor expansion coefficients at the knots,
# return a matrix (m+2) x (k+1) S filled with entries of the derivatives at the knots.
frlc=function(FR,LC,Vd){
  k=length(FR)-1
  m=length(LC)-2
  if(FR[k+1]!=LC[1])
  {
    stop("The last entry of the 'first row' argument is not the same as the first entry of the 'last column' argument.")
  }

  S=matrix(0,ncol=k+1,nrow=m+2)
  S[1,]=FR
  S[,(k+1)]=LC

  #To assure that illogical treatment of vectors and matices in R is fixed.
  if(min(dim(as.matrix(Vd)))==1){Vd=matrix(Vd,nrow=1)}
  ###############
  if(dim(Vd)[1]==1){
    if((k+1)!=dim(Vd)[2]){stop("Incompatible dimensions of the input matrices.")}
    for(i in 2:(m+2)){
      PS=S[i-1,]%*%toeplitz_lower(Vd)
      S[i,1:k]=PS[1:k]
    }
  }else{
    for(i in 2:(m+2)){
      PS=S[i-1,]%*%toeplitz_lower(Vd[(i-1),])
      S[i,1:k]=PS[1:k]
    }
  }
  return(S)
}

#First row -- last row verification as described in the paper
# FR -- (k+1) numeric vector. The first row in the output matrix
# LR -- (k+1) numeric vector. The last column in the output matrix
# Vd -- a (m+1) x (k+1) matrix if non-equidistant knots, or 1 x (k+1) if the equidistant knots (vectors are fine)
#      columnwise vectors  of the Taylor expansion coefficients at the knots,
# neqknots -- optional parameter utilized in the equidistant knots case for the number of internal knots,
#          default is zero, i.e. no internal knots
#
# return a matrix (m+2) x (k+1) U filled with entries of the derivatives at the knots.
# m should be less or equal k, if m>k, then the first k+1 rows of Vd are taken and the rest is desregarded.
# The message is generated and the output is the unique and not over-specified (k+2) x (k+1) solution U.
# If Vd 1 x (k+1), then m in the program is set to k to obtain the unique and not over-specified solution with k+2 knots.
# If 1 < m < k, the first k-m entries in LR, i.e. LR[1:(k-m)] are changed to guarantee a proper solution,
# the message is generated.
frlr=function(FR,LR,Vd,neqknots=0){

  if(length(FR)!=length(LR)){stop("Incompatible dimensions of the input matrices.")}

  k=length(FR)-1

  #To assure that illogical treatment of vectors and matices in R is fixed.
  if(min(dim(as.matrix(Vd)))==1){
    Vd=matrix(Vd,nrow=1)
    m=neqknots
  }
  ###############
  if(dim(Vd)[2]!=k+1){stop("Incompatible dimensions of the input matrices.")}


  if(dim(Vd)[1]>1){
    m=dim(Vd)[1]-1
    if(m>k){
      cat("\nThere are", m+2 ,"knots only the first", k+2, "knots are considered.\n\n")
      m=k
      Vd=Vd[1:(m+2),]
    }else{
      if(m<k){
        cat("\nThere are less than", k+2, "knots, the first",k-m,"entries of the", k+2,"nd row counting from the end in the input will be changed in the output.\n\n")
      }
    }
  }


  U=matrix(0,ncol=k+1,nrow=m+2)
  U[1,]=FR
  U[m+2,]=LR
  if(m>0){
    if(dim(Vd)[1]==1){
      A=toeplitz_lower(Vd)
      c=A[k+1,(k-m+1):k]
      A=A[(k-m+1):k,(k-m+1):k]
      D=rbind(A,c)
      E=matrix(0,ncol=m,nrow=m)
      for(i in 1:m){
        D=D%*%A
        E[m+1-i,]=c
        c=c%*%A
      }
    }else{
      A=toeplitz_lower(Vd[1,])
      D=A[(k-m+1):(k+1),(k-m+1):k]
      E=matrix(0,ncol=m,nrow=m)
      AA=diag(m)
      for(r in (m+1):2){
        A=toeplitz_lower(Vd[r,])
        c=A[k+1,(k-m+1):k]
        E[r-1,]=c%*%AA
        AA=A[(k-m+1):k,(k-m+1):k]%*%AA
      }
      D=D%*%AA
    }
    E=solve(E)
    U[2:(m+1),k+1]=(LR[(k-m+1):k]-FR[(k-m+1):(k+1)]%*%D)%*%E
  }

  U=frlc(U[1,],U[,(k+1)],Vd)
  return(U)
}

#----------------------------------------------------------------------------------------#

numToWord = function(x){
  if(x==1){
    return("First")
  } else if(x==2){
    return("Second")
  } else if(x==3){
    return("Third")
  } else{
    return(paste0(x, 'th'))
  }
}
#####
typeWords = function(t, w = TRUE){
  type = c("bs", "gsob", "twob", "spnt", "ndspnt", "sp")
  words = c(
    "Bspline functions\n", "gram-schmidt orthognal basis\n", "twosided orthognal basis\n",
    "dyadic splinet\n", "not fully dyadic splinet\n", "splines\n"
  )
  size = c(
    "basis functions\n", rep("orthogonal splines\n", 4), "spline functions\n"
  )
  if(w){
    res = words[which(type == t)]
  } else{
    res = size[which(type == t)]
  }
  return(res)
}

####
print.class <- function(object){
  # 1) type
  if(is.null(object@type)){
    cat("Splinets object\n")
  } else{
    cat(typeWords(object@type))
  }
  # 2) knots
  if(object@equid){
    cat(paste("Knots:", length(object@knots), "equaly distributed knots between",
              range(object@knots)[1], "and", range(object@knots)[2], "\n"))
  } else{
    cat(paste("Knots:", length(object@knots), "non-equaly distributed knots between",
              range(object@knots)[1], "and", range(object@knots)[2], "\n"))
  }
  # 3) size
  cat(paste("Size:", length(object@der), typeWords(object@type, w = FALSE)))
  # 4) order
  cat(paste("Order:", object@smorder, "\n"))
  invisible(object)
  # 5) support
  cat(paste("Support:", supp_info(object@supp)))
}
################
# the function for extracting supports information for online printing object
supp_info = function(supp){
  d = length(supp)
  if(d == 0){
    supp_message = "The full support range for each spline."
  } else if(d > 0){
    res = numeric(d)
    for(i in 1:d){
      res[i] = dim(supp[[i]])[1]
    }
    res = unique(res)
    if(min(res) == 1 & max(res) == 1){
      supp_message = "Not the full range support, a single support interval for each spline."
    } else{
      supp_message = paste("Not the full range support, the number of support intervals varies from", min(res), "to", max(res))
    }
  } 
  return(supp_message)
}


#The following function for a given set of knots and the order of smootheness evaluates the Taylor coefficients
#usefull for evaluation of the values of the splines in-between knots and stores in (n+1) x (k+1) matrix as the output

taylor_coeff=function(knots,k){
  n=length(knots)-2
  delta=diff(knots)
  Ad=matrix(0,nrow=n+1,ncol=k+1)
  Ad[,1]=rep(1,n+1)
  if(k>0){
    for(j in 1:k){
    Ad[,j+1]=(Ad[,j]*delta)/j
    }
  }
  return(Ad)
}

#For a vector of values in 'Vd' the following function produces the lower Toeplitz matrix.
#The function works on the vector, column and row formats of 'Vd'
toeplitz_lower=function(Vd){
  k=length(Vd)-1
  A=matrix(0,ncol=(k+1),nrow=(k+1))

    for(j in 1:(k+1)){
      A[j:(k+1),j]=Vd[1:(k-j+2)]
    }

  return(A)
}

#For two vectors V1, V2 function applies a convolution style operation and returns a new vector V
#such that 'toeplitz_lower(V)=toeplitz_lower(V1)%*%toeplitz_lower(V2)

toeplitz_conv=function(V1,V2){
  k=length(V1)-1
  V=convolve(V1,rev(V2),type="o")
  return(V[1:(k+1)])
}

#----------------------------------------------------------------------------------------#

#---------------------------------#
#--- Sec 2 AuxFun for evspline ---#
#---------------------------------#

#Function for evaluate splines. It is the aux function for function 'evspline'
#object: Splinets object
#sID: integer, specify which spline will be evaluated.
#x: arguments at which the splines are 
#return a vector of spline values with the same length as t.
evaluate_spline <- function(xi,supp,k,S, x){
  n = length(xi)-2
  if(length(supp) == 0){#if supp is NULL (the length is zero), then the full support is assumed.
    supp = matrix(c(1,n+2), ncol = 2)
  } 
  # Step 1) transform der matrix to the onesided form
  S = sym2one(S, supp)
  # Step 2) preparing loop over supp-intervals
  d = dim(supp)[1]
  evv = numeric(length(x)) # the vector of evaluations
  dd = supp[,2]-supp[,1] # for extracting sub der matrix given intervals
  st = 1 # start row of der
  # Step 3) loop over sub intervals
  for(j in 1:d){
    # Step 3.1) extract the sub der matrix into SS
    supp_temp = supp[j, ]
    en = st + dd[j]
    SS = S[st:en, ,drop=FALSE]
    st = en+1
    # Step 3.2) extract the corresponding part of the knots under supp_temp
    xi_supp = xi[supp_temp[1]:supp_temp[2]]
    # Step 3.3) extract the arguments within supp_temp
    id_l=head(which(x>=xi[supp_temp[1]]),1) # id of the first value in t within supp_temp
    if(length(id_l)!=0){ #if there are any points from the input in the particual interval
      id_r=tail(which(x<=xi[supp_temp[2]]), 1) # id of the last value in t within supp_temp
      if(length(id_r)!=0){ #if there are any points from the input in the particual interval
       if(id_l<=id_r){  #if there are any points from the input in the particual interval
         x_temp = x[id_l:id_r] # arguments within supp_temp
         # Step 3.4) find the intervals (xi_i, xi_{i+1}) that arguments x_temp fall in
         id_interval = findInterval(x_temp, xi_supp, all.inside = TRUE)
         # Step 3.5) calculate x-xi
         re_t = x_temp-xi_supp[id_interval]
         # Step 3.6) calculate (x-xi)^k/k! and save in T_matrix
         T_matrix = matrix(1, nrow=length(x_temp), ncol = k+1)
         if(k>0){
           for(i in 1:k){
           T_matrix[,(i+1)] = re_t*(T_matrix[,i])/i
           } 
         } 
         # Step 3.7) calculate evaluations, i.e. sum_{i=0}^k f^{(i)}(x-xi)^i/i!
         evv[id_l:id_r] = (T_matrix*SS[id_interval,])%*%rep(1,k+1)
      }
    }

      
    } #if there are any points from the input in the particual interval
  }
  
  return(evv)
}

#--------------------------------#
#--- Sec 3 AuxFun for lincomb ---#
#--------------------------------#

# Function for calculating the sum of two splines
# S1, S2: der matrix with one side form
# supp1, supp2: supports
# return the der matrix and supp of the sum of two splines
sum_two = function(S1, S2, supp1, supp2){
  supp = supp_union_bin(supp1, supp2) #the support of the union
  id=NULL
  for(i in 1:dim(supp)[1]){
    id=c(id,supp[i,1]:supp[i,2])
  } #the indexes in the support
  S = matrix(0,ncol=dim(S1)[2],nrow=length(id))
  id1=NULL
  for(i in 1:dim(supp1)[1]){
    id1=c(id1,supp1[i,1]:supp1[i,2])
  }
  id2=NULL
  for(i in 1:dim(supp2)[1]){
    id2=c(id2,supp2[i,1]:supp2[i,2])
  }
  # get the id of knots within the intersection and set differences.
  i_12 = intersect(id1, id2) # supp1 n supp2
  i_1 = setdiff(id1, i_12) # supp1\supp2
  i_2 = setdiff(id2, i_12) # supp2\supp1
  if(length(i_12) != 0){# der S for supp1 n supp2
    S[which(id %in% i_12),] = S1[which(id1 %in% i_12), ] + S2[which(id2 %in% i_12), ]
  }
  if(length(i_1 != 0)){# der S for supp1\supp2
    S[which(id %in% i_1),] = S1[which(id1 %in% i_1), ] 
  }
  if(length(i_2) != 0){# der S for supp2\supp1
    S[which(id %in% i_2),] = S2[which(id2 %in% i_2), ] 
  }
  return(list(S = S, supp = supp))
}

# function for doing linear combination iteratively by 'sum_two'
# S: list of der matrix
# supp: list of support matrix
# A: transformation matrix
# n: the number of inner knots
lincomb_supp = function(S, supp, A, n){
  n_l = dim(A)[1] # total number of linear combinations
  rS = list() # the resulting list of der matrix
  rsupp = list() # resulting matrix of support
  
  d = length(S)
  full_supp = matrix(c(1,n+2), ncol = 2) # generate the full range support
  if(length(supp) == 0){ # useful when all splines are full support, i.e. supp is empty set.
    supp = rep(list(full_supp), d) # this creates the list with all the elements the same
  }
  # switch to onesided representation
  for(i in 1:length(S)){
    S[[i]] = sym2one(S[[i]], supp[[i]])
  }
  # do linear combination one spline by one spline-
  for(i in 1:n_l){
    a = A[i, ]
    id = which(a != 0) # id of included splines
    n_i = length(id) # total number of included splines
    
    rS[[i]] = a[id[1]]*S[[id[1]]]
    rsupp[[i]] = supp[[id[1]]]
    
    if(n_i > 1){
      for(j in 2:n_i){
        temp = sum_two(S1 = rS[[i]], S2 = a[id[j]]*S[[id[j]]], supp1 = rsupp[[i]], supp2 = supp[[ id[j] ]])
        rS[[i]] = temp$S
        rsupp[[i]] = temp$supp
      }
    }
  }
  # switch back to symmetric representation
  for(i in 1:n_l){
    rS[[i]] = sym2one(rS[[i]], rsupp[[i]], inv = TRUE)
  }
  
  # fix support
  all_full = TRUE
  ct = 1
  while(all_full & ct <= length(rsupp)){
    all_full = is.full(rsupp[[ct]], full_supp)
    ct = ct + 1
  }
  if(all_full){
    rsupp = list()
  }
  
  return(list(S = rS, supp = rsupp))
  
}

# Function for linear combination with full der matrix
# S: a list of der matrix with sym representation
# supp: matrix of support, 2 columns
# k: order
# n: number of inner knots
# A: matrix
# return: 
#   1) a list of der matrices of linear combinations
#   2) a list of corresponding support matrices
lincomb_full = function(S, supp, k, n, A, SuppCont){
  l = dim(A)[1] # the number of combinations
  d = length(S) # the number of functions
  # fix the case with empty support case
  full_supp = matrix(c(1,n+2), ncol = 2) # generate the full range support
  if(length(supp) == 0){ # useful when all splines are full support
    supp = rep(list(full_supp), d)
  }
  #------------------------------------------------------------------------------------------------#
  # Step 1: stack all full der matrix in a matrix for calculating the linear combination
  temp_matrix = matrix(0, nrow = (n+2)*(k+1) , ncol = d) 
  for(i in 1:d){
    # 1.1: transform from sym 2 onesided representation
    temp_S = sym2one(S[[i]], supp[[i]])
    # 1.2: retrive the full der matrix if the support is not full support
    if(!is.full(supp[[i]], full_supp)){
      temp_S = full_der(temp_S, supp = supp[[i]], n = n) # get the full der matrix
    }
    # 1.3: stack into matrix
    temp_matrix[,i] = c(temp_S) # get the one side representation
  }
  # Step 2: linear combination by big matrix calculation
  temp_matrix = temp_matrix %*% t(A)
  # Step 3: retrive results
  S_new = supp_new = list()
  if(SuppCont){
    for(i in 1:l){
      S_new[[i]] = sym2one(matrix(temp_matrix[,i], byrow = FALSE, nrow = n+2, ncol = k+1), inv = TRUE)
      temp_res = exsupp(S_new[[i]])
      S_new[[i]] = temp_res$rS
      supp_new[[i]] = temp_res$rsupp
    }
  } else{
    for(i in 1:l){
      S_new[[i]] = sym2one(matrix(temp_matrix[,i], byrow = FALSE, nrow = n+2, ncol = k+1), inv = TRUE)
    }
  }
  
  #------------------------------------------------------------------------------------------------#
  # fix support, if supp is full then set it as empty
  all_full = TRUE
  ct = 1
  while(all_full & ct <= length(supp_new)){
    all_full = is.full(supp_new[[ct]], full_supp)
    ct = ct + 1
  }
  if(all_full){
    supp_new = list()
  }
  #------------------------------------------------------------------------------------------------#
  return(list(S=S_new, supp=supp_new))
}

#---------------------------------#
#--- Sec 4 AuxFun for integral ---#
#---------------------------------#

# knots: knots vector 
# supp: a vector, support of input S
# taylor: taylor coefficients matrix
# S: derivative matrix
dint_single = function(knots, supp, taylor, S){
  l1 = supp[1]
  l2 = supp[2]
  S = sym2one(S, supp)
  S = S[-dim(S)[1], ]
  
  k = dim(S)[2]-1
  taylor = cbind(taylor, diff(knots)^(k+1)/factorial(k+1))
  taylor = taylor[l1:(l2-1), -1]
  res = sum(apply(taylor*S,1,sum))
  return(res)
}

dint_engine = function(knots, supp, taylor, S){
  d =  dim(supp)[1]
  SS = der_split(S, supp)
  res = dint_single(knots, supp[1,,drop=FALSE],taylor,SS[[1]])
  if(d > 1){
    for(i in 2:d){
      res = res + dint_single(knots, supp[i,,drop=FALSE],taylor,SS[[i]])
    }
  }
  return(res)
}

# Function for calculate indefinite integral over a single component 
# taylor: (n+1) x (k+2) taylor coefficients matrix for the set of knots corresponding to S. 
# It has one more column than 'S' for the integral derivatives
# 
# S: (n+2) x (k+1) derivative matrix with sym representation made of a single component having n+2 knots
# c_0: the value at the left-hand side end
# output: (n+2) x (k+1) derivative matrix for the indefinite integral in the sym representation
# it has to be noted that the first column does not have typically value zero in the bottom position


int_single = function( taylor, S, c_0 = 0){
  
  S = sym2one(S) # transform to one sided it is expected that 'S' corresponds to a single component 
  n = dim(S)[1]-2
  S = cbind(numeric(n+2), S) # here all the values of the derivatives of the integral are set 
  S[1,1] = c_0 #the value of the integral at the first knot
  # start to fix the first column, i.e. values of the integral
  
  for(j in 2:(n+2)){
    S[j,1] = sum(S[j-1,]*taylor[(j-1), ])
  }
  
  S = sym2one(S,inv =TRUE) #S[(n+2),1] typically is non-zero
  return(S)
}


#This function computes the coefficients for the matrix of derivative of the integral of a spline with 
#arbitrary number of the components in the support
#supp - Nsupp x 2 matrix of support indexes usually taken from a single spline in the Splinet object 
#       NULL support is not covered by the function, it assumes the full range in this case.
#       the values in supp are from 1:(n+2)
#taylor -(n+1) x (k+2) taylor expansion coefficients at the knots for the integral result 
#S the  matrix of derivatives values in the symmetric representation and corresponding to supp thus the 
#  total size is m x (k+1) where m=sum(supp[,2]-supp[,1]+1)
#returns the corresponding matrix of the derivatives m x (k+2) for the integral in the symmetric representation
# on the componnents of the support. The boundary conditions are not satisfied but are the correct values of derivatives. 


int_engine = function(supp, taylor, S ,epsilon = 1e-7){
   #epsilon -- the parameter to determine if one deals with the spline after integration to
   #extract the support using exsupp
  Nsupp = dim(supp)[1] #The number of elements in the support (null case is not possible )
  
  rS = NULL #space for the computed der matrix of the integral
 
  
  #Evaluation of the matrix of the derivatives over the  components of the support
  c_0=0
  beg=1 #evaluation of the current beginning in S corresponding to the component (should be equal to one)

    for(i in 1:Nsupp){
      
      l1 = supp[i,1]
      l2 = supp[i,2]
      len=l2-l1
      
      loc_taylor = taylor[l1:(l2-1),] #extracting part of the taylor relevant for the support
      
      AA = int_single(loc_taylor, S[beg:(beg+len),,drop=FALSE],c_0=c_0)
      
      rS=rbind(rS,AA) #updating the output matrix 
      
      c_0=rS[beg+len,1] #The RHS is the value of the spline at the last knot in the component
                        #It should be the balue for the value of the integral at the beginning of 
                        #the next component
      beg=beg+len+1 #the beginning for the next loop
      
    }
  
  return(rS)
}

#--------------------------------#
#--- Sec 5 AuxFun for gramian ---#
#--------------------------------#

# function for calcalting the innder product between two splines with single interval support
# der1/2: der matrix
# supp1/2: support, vectors
# k: order
# n: number of inner knots

inner_single <- function(der1, der2, supp1, supp2, FF, D, n, k){
  supp_1 = supp1[1]:supp1[2]
  supp_2 = supp2[1]:supp2[2]
  supp = setdiff(intersect(supp_1, supp_2), n+2)
  l = length(supp)
  if(l > 0){
    S1 = der1[supp-supp1[1]+1, ,drop=FALSE]
    S2 = der2[supp-supp2[1]+1, ,drop=FALSE]
    d = l-1
    Conv_S = numeric(2*k+1)
    FF_temp = FF[supp, ,drop=FALSE ]
    for(i in 1:(d+1)){
      A=FF_temp[i,]*S1[i,]
      B=FF_temp[i,]*S2[i,]
      B=B[(k+1):1]
      Conv_S = Conv_S + convolve(A, B, type = "open")
    }
    res = Conv_S%*%D
  } else{ res = 0 }
  return(res)
}


# function for calculating the innder product between two splines
# der1, der2: der matrix
# supp1, supp2: Support matrix
inner_engine = function(der1, der2, supp1, supp2, FF, D, n, k){
  d1 = dim(supp1)[1]
  d2 = dim(supp2)[1]
  S1 = der_split(der1, supp1)
  S2 = der_split(der2, supp2)
  res = 0
  for(i in 1:d1){
    for(j in 1:d2){
      ss = intersect(supp1[i,1]:supp1[i,2], supp2[j,1]:supp2[j,2])
      if(length(ss) > 1){
        res = res + inner_single(S1[[i]], S2[[j]], supp1[i, ], supp2[j, ], FF, D, n, k)
      }
    }
  }
  return(res)
}
#The following function evaluates an n+1 x k+1 matrix
#FF_{r,j} dxi_r^(j-1/2)/(j-1)!, r=1,...,n+1; j=1,...,k+1, this matrix appears in the arkiv paper
#on splinets, but in the other R-paper on the arkiv there is a mistake in Proposition 3. 
FF_Matrix <- function(knots, k){
  n = length(knots) - 2
  FF=matrix(0,ncol=k+1,nrow=n+1)
  df_xi=as.matrix(diff(knots))
  FF[,1]=sqrt(df_xi) #see the formula on the splinets
  for(j in 1:k){FF[,j+1]=FF[,j]*df_xi/j}
  return(FF)
}

# function for efficient calculating the gram matrix for bspline basis functions
# The inputs must be the components of a set of bspline basia functions 
# IMPORTANT: It is assumed that the B-splines are normalized or wrong results are returned!!!
# xi, knots
# k, order
# S, list of der matrices
# supp, a list of support
# return: gram matrix of so, which is a band matrix
bandmatrix = function(xi, k, S, supp){
  n_so = length(S)
  H = diag(n_so)
  n = length(xi) - 2
  FF = FF_Matrix(xi, k) #according to the formula on the inner product in the original arkiv paper
  D = as.matrix(1/(1:(2*k+1)))
  for(i in 1:n_so){
    S[[i]] = sym2one(S[[i]], supp[[i]])
  }
  if(n_so-k>0){# it may happen that n_so-k=0: we have on one group on one level
  for(i in 1:(n_so-k)){
    for(j in (i+1):(i+k)){
     # H[i,j] = H[j,i] = inner_engine(der1 = S[[i]], der2 = S[[j]],
     #                                supp1 = supp[[i]], supp2 = supp[[j]],
     #                               FF = FF, D = D, n = n_so+k-1, k = k)
       H[i,j] = H[j,i] = inner_engine(der1 = S[[i]], der2 = S[[j]],
                               supp1 = supp[[i]], supp2 = supp[[j]],
                               FF = FF, D = D, n = n, k = k)
    }
  }
}
  if(k != 1){
    for(i in (n_so-k+1):(n_so-1)){
      for(j in (i+1):n_so){
        #H[i,j] = H[j,i] = inner_engine(der1 = S[[i]], der2 = S[[j]],
        #                               supp1 = supp[[i]], supp2 = supp[[j]],
        #                               FF = FF, D = D, n = n_so+k-1, k = k)
        H[i,j] = H[j,i] = inner_engine(der1 = S[[i]], der2 = S[[j]],
                                        supp1 = supp[[i]], supp2 = supp[[j]],
                                       FF = FF, D = D, n = n, k = k)
      }
    }
  }
  return(H)
}

#----------------------------------#
#--- Sec 6 AuxFun for 'bspline' ---#
#----------------------------------#

#The following four functions are created for generate the der matrix
lambda1 <- function(xi, k){
  n = length(xi) - 2
  l = (n-k+1)*(k+1)
  res = array(0, dim = c(l, 1))
  
  xi = head(xi, -1)
  
  for(ord in 1:k){
    vec = seq(ord+1, l, by = k+1)
    res[vec] = diff(xi, ord)[1:length(vec)]
  }
  return(res)
}

lambda2 <- function(xi, k){
  n = length(xi) - 2
  l = (n-k+1)*(k+1)
  res = array(0, dim = c(l, 1))
  
  for(ord in 1:(k+1)){
    vec = seq(ord, l, by = k+1)
    res[vec] = tail(-diff(xi, k+2-ord), length(vec))
  }
  return(res)
}

coeff1 <- function(xi, k){
  matrix(1/c(matrix(rep(diff(head(xi, -1), k), k+1), byrow = TRUE, nrow = k+1)), ncol = 1)
}

coeff2 <- function(xi, k){
  matrix(1/c(matrix(rep(-diff(tail(xi, -1), k), k+1), byrow = TRUE, nrow = k+1)), ncol = 1)
}

#------------------------------#
#--- Sec 7 AuxFun for grsch ---#
#------------------------------#

# Implement the generic function 1 (algorithm 1): gram-schmidt orthonormalization
# It is worth to remember that if the nb x nb output P is orthogonalization in the problem
# where nb x nb input A=I, then k x nb output AP is the output if the input is k x nb A, or, 
# grscho(A,H)=A%*%grscho(I,H) 
# # Little example:
# I=diag(20)
# X=matrix(rnorm(20*20),nrow=20)
# H=t(X)%*%X
# A=matrix(rnorm(40*20),nrow=40)
# 
# P=grscho(I,H)
# PA=grscho(A,H)
# 
# PA-A%*%P


grscho=function(A, H){
  nb=dim(H)[1]
  B=A #the k x nb, k>=nb matrix output with columns being orthonormalized columns of A
  B[,1]=A[,1, drop=FALSE]/sqrt(H[1,1]) #the first output vector normalized
  for(i in 1:(nb-1)){
    A[,(i+1):nb]=(A[,(i+1):nb]-(A[,i, drop=FALSE]%*%H[i,(i+1):nb,drop=FALSE])/H[i,i])
    H[(i+1):nb,(i+1):nb]=(H[(i+1):nb,(i+1):nb]-t(H[i,(i+1):nb,drop=FALSE])%*%H[i,(i+1):nb,drop=FALSE]/H[i,i])
    B[,i+1]=A[,i+1,drop=FALSE]/sqrt(H[i+1,i+1]) #subsequent column to the output
  }
  B
}

#---------------------------------#
#--- Sec 8 AuxFun for twosided ---#
#---------------------------------#

# Implement the generic function 2 (algorithm 2): Symmetrized orthonormalization
# x: 2 columns matrix contains coefficients w.r.t. the basis
# h: the inner product of two vectors w.r.t. the basis
symo = function(x, h){
  a1 = (1/sqrt(1+h) + 1/sqrt(1-h))/2
  a2 = (1/sqrt(1+h) - 1/sqrt(1-h))/2
  res = x
  res[,1] = a1*x[,1]+a2*x[,2]
  res[,2] = a2*x[,1]+a1*x[,2]
  return(res)
}

# Implement the generic function 3 (algorithm 3): Symmetrized Gram-Schmidt orthonormalization
# The same comments as in grscho apply and the following example illustrates it
# If A is NULL that the problem is solved for the identity I
# It exploits the following simple while not immediately obvious relation:
# If the nb x nb output P is orthogonalization in the problem
# where nb x nb input A=I, then k x nb output AP is the output if the input is k x nb A.

sgrscho = function(A = NULL , H){
  nb = dim(H)[1]; np = floor(nb/2) # the number of vectors and pairs
  
  A_input=A #used at the end to get a result for the case if A is not NULL
  
  A=diag(nb)
  B = A
  
  BR = matrix(0, ncol = 2*np, nrow = nb)
  HR = matrix(0, ncol = 2*np, nrow = 2*np)
  J = c(2*(1:np)-1, 2*(1:np))
  K = c(nb:(nb-np+1), 1:np)
  BR[, J] = A[, K]
  HR[J,J] = H[K,K]
  
  BL = matrix(0, ncol = nb, nrow = nb)
  HL = matrix(0, ncol = nb, nrow = nb)
  J = c(J, nb)
  K = c(1:np, nb:(nb-np+1), np+1)
  BL[, J] = A[, K]
  HL[J,J] = H[K,K]
  
  BL = grscho(BL, HL)
  BR = grscho(BR, HR)
  
  B[, np+1] = BL[, nb] # center for the odd case
  
  for(i in 1:np){
    X = cbind(BL[, 2*i-1], BR[, 2*i-1])
    h = (X[, 1]%*%H%*%X[, 2, drop=FALSE])[1,1]
    B[, c(i, nb-i+1)] = symo(X, h)
  }
  if(!is.null(A_input)){
    if(dim(A_input)[2]!=dim(H)[1]){cat("The dimension of the input matrices do not match.\n
                                 The computations assume that the first input in
                                 symmetric orthonormalization algorithm is NULL.")}else{
                                   B=A_input%*%B
                                   } 
                                 }
  return(B)
}

#--------------------------------#
#--- Sec 9 AuxFun for splinet ---#
#--------------------------------#

# The core recursive step for 'dyadicsplinet' function.
# It is implementation of the recursion discussed in the paper on splinets.
# A - d x d matrix, a columnwise representation of d linearly independent vectors in a certain basis made
# of d elements
# H - d x d matrix, Gram matrix for the basis in which vectors in A are represented in 
#     its columns; It is a band matrix with the band width equal to 2k-1.
# k - the band width parameter
# Toep - logical, to indicate if efficiency of the algorithm due to Toeplitz input matrix 
# should be used
# Depth - will be used for the non-dyadic Toeplitz (not implemented yet)
# The fully dyadic structure is assumed so that d=k(2^N-1), for some positive integer N

dyadic_engine = function(A, H, index, k, Toep = FALSE ){ 
  d = length(index)
  id_grou = matrix(1:d, nrow = k) #The problem has n_grou k-tuples
  #each row corresponds to a k-tuple
  n_grou = d/k
  old_H = H
  # SGSO step -- the symmetric orthogonalization of the lowest level
  if(Toep==FALSE){
    if(k>1){
      # sym orth
      
      for(grou in seq(1, n_grou, by = 2)){ # grou = tuplet/group, by = 2, because 
        # every other k-tuple belongs to the lowest level
        temp_grou = id_grou[, grou]  # k dimensional vector with index of the k-tuple 'grou'
        spli = index[temp_grou]      # indexes in the original A matrix for the vectors in  the k-tuple 'grou'  
        temp_A = A[spli, spli]
        temp_H = H[spli, spli]
        
        A[spli, spli] = sgrscho(temp_A, temp_H) #GSO orthogonalization of the bottom level k-tuplets
      }
      # update Gram matrix
      for(grou in seq(2, n_grou, by = 2)){
        spli_C = index[id_grou[, grou]]
        spli_L = index[id_grou[, grou-1]]
        spli_R = index[id_grou[, grou+1]]
        H[spli_L, spli_C] = t(A[, spli_L])%*%H%*%A[, spli_C]
        H[spli_C, spli_L] = t(H[spli_L, spli_C])
        H[spli_R, spli_C] = t(A[, spli_R])%*%H%*%A[, spli_C]
        H[spli_C, spli_R] = t(H[spli_R, spli_C])
      }
    }
    # GSO step
    for(grou in seq(2, n_grou, by = 2)){ #Orthogonalizing the rest of the matrix with respect 
      #to the lowest level
      spli_C = index[id_grou[, grou]]
      spli_L = index[id_grou[, grou-1]]
      spli_R = index[id_grou[, grou+1]]
      
      H_L = H[spli_L, spli_C]
      H_R = H[spli_R, spli_C]
      A_L = -A[spli_L,spli_L]%*%(H_L)
      A_R = -A[spli_R,spli_R]%*%(H_R)
      
      A[spli_L, spli_C] = A_L
      A[spli_R, spli_C] = A_R
      
    }
  }else{
    if(k>1){ #This requires the original problem to by dyadic and Toeplitz
      #recall that if it is Toeplitz but not dyadic the matrix is complemented at 
      #the end to a complete dyadic case through diagonal blocks at the beginning
      #and at the end. This destroy the Toeplitz condition for the extended matrix.
      #In future, one can implement a different form of the orthogonalization for 
      #The non-dyadic case by spliting to a series of disjoing smaller dyadic problems and 
      #at the junctions of these problems applty the orthonormalization. This would
      #work well with the non-dyadic Toeplitz case. 
      
      #Sym orthogonalization for the first k-tuple, k>1
      temp_grou = id_grou[,1]  # k dimensional vector with index of the k-tuple 'grou'
      spli = index[temp_grou]      # indexes in the original A matrix for the vectors in  the k-tuple 'grou'  
      temp_A = A[spli, spli]
      temp_H = H[spli, spli]
      
      A[spli, spli] = sgrscho(temp_A, temp_H) #
      
      #All the remaining portion of the matrices corresponding to the lowest level will be the same
      if(n_grou>2){
        for(grou in seq(3, n_grou, by = 2)){ #These are indices of the blocks of A to which the above 
          #computed values has to be copied since due to the Toeplitz
          #property they will be identical
          temp_grou=as.vector(id_grou[,grou])        #k x (2^{N-1}-1) matrix of indices of the groups corresponding
          #to the block matrices in A to be modified,  
          #It has to be translated to a sequence by column in order
          #to provide proper indexes of the matrices
          spli2=index[temp_grou]
          A[spli2,spli2]=A[spli,spli] #repeating the first matrix at all other spots
        }
      }
      # update Gram matrix
      grou=2 #Furthest to the left 
      spli_C = index[id_grou[, grou]] #indexes for the matrices above the lowest level
      spli_L = index[id_grou[, grou-1]] #index for the LHS matrix at the lowest level
      spli_R = index[id_grou[, grou+1]] #index for the RHS matrix at the lowest level
      H[spli_L, spli_C] = t(A[, spli_L])%*%H%*%A[, spli_C] #updating of the inner products with respect the LHS vectors
      H[spli_C, spli_L] = t(H[spli_L, spli_C])
      H[spli_R, spli_C] = t(A[, spli_R])%*%H%*%A[, spli_C] #updating of the inner products with respect the RHS vectors
      H[spli_C, spli_R] = t(H[spli_R, spli_C])
      if(n_grou>3){
        for(grou in seq(4, n_grou, by = 2)){
          spli_C2 = index[id_grou[, grou]]
          spli_L2 = index[id_grou[, grou-1]]
          spli_R2 = index[id_grou[, grou+1]]
          H[spli_L2, spli_C2] = H[spli_L, spli_C]
          H[spli_C2, spli_L2] = H[spli_C, spli_L]
          H[spli_R2, spli_C2] = H[spli_R, spli_C]
          H[spli_C2, spli_R2] = H[spli_C, spli_R]
        }
      }
    }
    # GSO step 
    #Orthogonalizing the rest of the matrix with respect 
    #to the lowest level
    grou=2
    spli_C = index[id_grou[, grou]]
    spli_L = index[id_grou[, grou-1]]
    spli_R = index[id_grou[, grou+1]]
    
    H_L = H[spli_L, spli_C]
    H_R = H[spli_R, spli_C]
    A_L = -A[spli_L,spli_L]%*%(H_L)
    A_R = -A[spli_R,spli_R]%*%(H_R)
    
    A[spli_L, spli_C] = A_L
    A[spli_R, spli_C] = A_R
    if(n_grou>3){
      for(grou in seq(4, n_grou, by = 2)){ 
        spli_C2 = index[id_grou[, grou]]
        spli_L2 = index[id_grou[, grou-1]]
        spli_R2 = index[id_grou[, grou+1]]
        
        A[spli_L2, spli_C2] = A_L
        A[spli_R2, spli_C2] = A_R
        
      }
    }
  }
  
  
  H = diag(dim(A)[1])
  temp_id = index[c(id_grou[,seq(2, n_grou, by = 2)])]
  H[temp_id,temp_id] = t(A)[temp_id, ]%*%old_H%*%A[, temp_id]
  
  return(list(A = A, H = H))
}

# This is the main dyadic algorithm function that requires dyadic shape of the input
# 
# par H matrix, a dyadic non-negative band matrix to be orthogonalized
# par k integer, a width of the band
# par Toep logical, a flag to indicate when the Toeplitz matrix (constants along the diagonals)
#         are involved and thus faster algorithm can be used, in the present implementation t
#         the efficiency is implemented only if the original problem is the fully dyadic case. 
# return matrix P such that P^THP=I (or more general P^TA^THAP), it may return zeros at the end of the diagonal if 
#        the input was singular, but robustness of this has not been tested so non-singular
#        matrices are recommended for the input. In other words, columns AP yields orthogonalized
#        vectors in the basis in which columns of A represented the original vectors (A=I is often the case)
# the main algorithm, it effectively works on matrices and it utilizes dyadic engine
# The dyadic structure means that:
# A is d x d matrix, the columns of which represent some linearly independent vectors, 
# H is d x d k-band matrix of the inner vectors given in A
# d=k*(2^N-1), for some positive integer N, 
# The vectors in A (and thus H) are forming the following structure:
# They are put in 2^N-1 k-tuples. The bottom, 1st layer is made of every second k-tuple removed
# starting from the second and ending with the second last, so that there is the total 2^{N-1} k-tuples
# in this layers and 2^{N-1}-1 other k-tuples that can be put similarily in other layers
# l=2,...,N-1, each of them having 2^{N-l} k-tuples in it. 

dyadicsplinet  = function(A, H, k, Toep = FALSE ){
  d = dim(A)[1]
  
  #########
  #The regular check of the dyadic form
  N=log2((d+k)/k) #The number of layers plus 'one'
  if(N-floor(N)!=0) stop("The size of the input matrices does not have ther required dyadic structure.")
  #End of the check
  #########
  
  index = 1:d #index controls how many time the recursion is performed, if index==k,
  #then only the top layer made of one k-tuple is left
  
  while (length(index)>k){ #The engine has the whole matrices with all levels instead
    #only the ones that need to be modify (the `top` of the pyramid),
    #the argument `index` controls what will be modified inside the engine.
    #but it would be more natural to have this outside of the engine).
    #This structure could be modified in the future, so that
    #the output from the engine would be properly modified the 
    #input columns in dxd' A and the d'xd' H, where d=k*(2^N-1), 
    #d'=k*(2^(N-l+1)-1), where the current level of the orignal matrix processed
    #Technically speaking d can be of any form for the engine to work
    #   res = dyadic_engine2(A = A[,index], H = H[index,index],  k = k, Toep)  
    #Since as it stands now the structure of the code should not affect
    #the performance we will stick to having the entire H entering in the input 
    #and diagonal A (in a sense we do not use the option to put arbitrary A) 
    #performing computation to get P=res$A[index,index] such that for A=A[,index],
    #AP is the modified and basis with the bottom level at a given step orthogonalized. 
    # It exploits the following simple while not immediately obvious relation:
    # If the nb x nb output P is orthogonalization in the problem
    # where nb x nb input A=I, then k x nb output AP is the output if the input is k x nb A.
    
    
    
    #   A[, index] = (A[,index] %*% res$A) # AP transformed column vectors corresponding to index
    # for the lowest level it is ready to go to the output
    
    res = dyadic_engine(A = diag(d), H = H, index = index, k = k, Toep)
    #Returns res$A and res$H so that the vectors in A corresponding to the lowest layers
    #are orthonormalized, and the rest of the vectors are orthogonalized to the
    #one in the lowest layer and in such a way that their Gram matrix is still a band
    #matrix, while the Gram matrix corresponding to the lowest layer is the identity
    #these Gram matrices are kept in the output res$H 
    A[, index] = (A[,index] %*% res$A[index, index]) # AP transformed column vectors corresponding to index
    # for the lowest level it is ready to go to the output
    
    H = res$H
    id_grou = matrix(1:length(index), nrow = k)
    index = index[c(id_grou[, seq(2, dim(id_grou)[2], by = 2)])] #Reducing the size of the matrices 
    #through reducing the size 'index' by k2^N' from k(2^N'-1) 
    #to k(2^(N'-1)-1) 
  }
  if(k>1){ #This is symmetric orthogonalization of the top level made of one k-tuple
    spli = index #At this stage index is made of a single k-tuple
    temp_A = sgrscho(H=H[spli, spli]) #The default for the input A is NULL and then diagonal is used inside 'sgrscho'
    temp = diag(d)
    temp[spli,spli] = temp_A
    A[, spli] = A[, spli]%*%temp[spli, spli] #modifying the k-entries in A
  }
  return(A)
}



# This is the main dyadic algorith function 
# do the dyadic orthogonalization algorithm -- the main algorithm, it works  on matrices
# and utilizes dyadic engine but does not require dyadic shape of the input
# par H matrix, a non-negative band matrix to be orthogonalized
# par k integer, a width of the band
# par Toep logical, a flag to indicate when the Toeplitz matrix (constants along the diagonals)
#         are involved and thus faster algorithm can be used
# return matrix P such that P^THP=I, it may return zeros at the end of the diagonal if 
#        the input was singular, but robustness of this has not been tested so non-singular
#        matrices are recommended for the input. 


dyadiag = function(H, k, Toep){
  #computing gram matrix for the input Bsplines and padding 
  #with with 'ones' on the both end of the diagonal to make it fully dyadic
  n=dim(H)[1] + k - 1
  N=log2((n+1)/k)
  if(N-floor(N)!=0){#the non-dyadic case
    Toep=FALSE #The Toeplitz efficiency not implemented for non-dyadic case. The main reason
    #is difficulty of keeping track in the dyadic algorithm which terms should be adjusted
    #do not need to be evaluated. If the `Toeplitz efficiency` is to be implemented for non-dyadic case
    #then another algorithm should be used: spliting the main case into the decreasing maximal dyadic
    #cases and applying the full dyadic algorithm with Toeplitz efficiency to each of them and then 
    #dealing with the k-tuples in between the dyadic components. Work for future. 
    N=ceiling(N)  #number of the levels
    n_aug=k*(2^N-1) #the size of the extended matrix: 'k*2^N-1' is the number internal knots,
    # k*2^N-1-k+1 number of basis elements = the dimension of augmented gram
    aug_H=diag(1,n_aug)
    n_U = floor((n_aug-dim(H)[1])/2) #The number of extra at the beginning
    n_D = n_aug-dim(H)[1]-n_U        #The number of extra at the end
    aug_H[(n_U+1):(n_aug-n_D), (n_U+1):(n_aug-n_D)]=H
    P = dyadicsplinet(A = diag(n_aug), H = aug_H, k = k,Toep) #this is the main algorithm
    #for dyadic orthogonalization see AuxFun.R for its code
    P = P[(n_U+1):(n_aug-n_D), (n_U+1):(n_aug-n_D)]
  }else{
    P = dyadicsplinet(A = diag(dim(H)[1]), H = H, k = k,Toep) #this is the main algorithm
    #for dyadic orthogonalization see AuxFun.R for its code 
  }
  return(P)
}


# function for calculating the dyadic number given the number of input bspline bases function.
# assist  'net_structure'.
aug_bound = function(n_so, k){
  # algo 4.4: given the number of splines and order, return the relative position of splines in the
  # augmented Gram matrix, that is (n_U+1):(n_so_1-n_D). And the number of support levels
  n = n_so + k + 1 - 2
  N = ceiling(log((n+1)/k)/log(2))
  n_U = floor((k*2^N-n-1)/2)
  n_D = k*2^N-n-1-n_U
  return(list(n_U = n_U, n_D = n_D, N = N))
}
# Generate the map between sequential order and net structure
# n_so: the number of splines
# k: smorder
# return: a matrix, from left to right, each column present 
#   1) 'Seq_ID': the sequential id
#   2) 'supp_level': the support level, the top of the dyadic structure corresponds to the smallest level (1)
#   3) 'tuplet_ID': the id of tuplet where the function stay (relatively to the fully dyadic structure).
#   4) 'basis_ID': the id of basis element within tuplet
net_structure = function(n_so, k){
  bound = aug_bound(n_so, k)
  temp = c(rep(0, bound$n_U), 1:n_so ,rep(0, bound$n_D))
  temp = matrix(temp, nrow = k)
  temp_res = array(0, c(dim(temp), 4))
  temp_res[,,1] = temp
  n_tuplets = dim(temp_res)[2]
  # supp level
  index = 1:n_tuplets
  supp_level = bound$N
  while(length(index) != 0 ){
    temp_id = index[seq(1,length(index), 2)]
    temp_res[,temp_id,2] = supp_level
    supp_level = supp_level-1
    index =  setdiff(index, temp_id)
  }
  # tuplet
  temp_res[,,3] = matrix(rep(1:n_tuplets, k), byrow = TRUE, nrow = k)
  # basis
  temp_res[,,4] = matrix(rep(1:k, n_tuplets), nrow = k)
  # transform
  res = apply(temp_res, 3, c)
  res = res[which(res[,1] != 0), ]
  colnames(res) = c("Seq_ID", "supp_level", "tuplet_ID", "basis_ID")
  return(res)
}

#------------------------------#
#--- Sec 10 Auxfun for plot ---#
#------------------------------#

plot.spline = function(object, x = NULL, sID = NULL, vknots=TRUE, mrgn=2, type='l', bty="n",col='deepskyblue4',
                       lty=1, lwd=2, xlim=NULL, ylim = NULL, xlab="", ylab = "",  ...){
  if(length(x) == 1){
    stop("x should not be a single value") #if x = NULL sampling the arguments will be made in evaluate function by the default sampling rate
  }
  if(length(sID) == 0){
    sID=1:length(object@der)
  }
  Nsp=length(sID)
  y = evspline(object, x = x, sID = sID)
  Arg=y[,1]
  Val=y[,-1,drop=FALSE]
  
  if(is.null(xlim)){
    xlim = range(Arg)
  }
  if(is.null(ylim)){
    ylim = range(Val)
  }
  par(mar = c(mrgn, 2*mrgn, mrgn, 2*mrgn))
  plot(Arg,Val[,1],type=type,bty=bty,col=col,xlim=xlim,ylim=ylim,
       xlab=xlab,ylab=ylab,lty=lty,lwd=lwd,...)
  if(Nsp>1){
    ourcol=c( 'darkorange3', 'goldenrod', 'darkorchid4',
             'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue','deepskyblue4')
    for(i in 2:Nsp){
      lines(Arg,Val[,i],col=ourcol[(i-2)%%8+1],lty=lty,lwd=lwd,...)
    }
  }
  if(vknots){
    abline(v = object@knots, lty = 3, lwd = 0.5)
  }
  abline(h = 0, lwd = 0.5)
}

######

lines.spline = function(object, x = NULL, sID = NULL, type='l', bty="n",col='deepskyblue4',
                        lty=1, lwd=2){
  if(length(x) == 1){
    stop("x should not be a single value")
  }
  if(length(sID) == 0){
    sID=1:length(object@der)
  }
  
  Nsp=length(sID)
  y = evspline(object, x, sID = sID)
  Arg=y[,1]
  Val=y[,-1, drop = FALSE]
  
  
  lines(Arg,Val[,1],col=col,lty=lty,lwd=lwd)
  if(Nsp>1){
    ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
             'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
    for(i in 2:Nsp){
      lines(Arg,Val[,i],col=ourcol[i%%8+1],lty=lty,lwd=lwd)
    }
  }
}
#########
plot.basis = function(object){
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  n_so = length(object@der)
  y = evspline(object)
  plot(y[,1], y[,2], type = "l", main = paste(numToWord(object@smorder), "order B-spline basis"),
       xlab = "", ylab = "", lwd = 2, col = ourcol[1], ylim = range(y), bty="n")
  if(n_so > 1){
    for(i in 3:(n_so+1)){
      points(y[,1], y[,i], type = "l", lwd = 2, col = ourcol[i%%8+1])
    }
  }
  abline(h = 0)
  xi = object@knots
  abline(v = xi, lty = 3, lwd = 0.5)
}
#############
plot.obasis = function(object){
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  n_so = length(object@der)
  y = evspline(object)
  plot(y[,1], y[,2], type = "l", main = paste(numToWord(object@smorder), "order", typeWords(object@type)),
       xlab = "", ylab = "", lwd = 2, col = ourcol[1], ylim = range(y), bty="n")
  for(i in 2:n_so){
    points(y[,1], y[, i+1], type = "l", lwd = 2, col = ourcol[i%%8+1])
  }
  xi = object@knots
  abline(v = xi, lty = 3, lwd = 0.5)
  abline(h = 0)
}
################
plot.splinet = function(object,lwd=2,mrgn=2){
  n_so = length(object@der)
  xi = object@knots
  k = object@smorder
  n = length(xi)-2
  y = evspline(object)
  # plot setups
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  # margin_height = 3
  par(mar = c(mrgn, 2*mrgn, mrgn, 2*mrgn))
  net_str = net_structure(n-k+1, k)
  n_level = max(net_str[, 2])
  layout(matrix(1:n_level, n_level, 1))
  
  
  for(i in 1:n_level){
    seqID = net_str[which(net_str[,2] == i), 1]
    
      plot(y[, 1], y[, seqID[1]+1], type = "l", ylab = "", bty = "n",
           ylim = range(y[, seqID+1]), xlab = "", col = ourcol[seqID[1]%%8+1], lwd = lwd,
           main = paste(numToWord(k), "order splinet basis" )) 
    
    for(j in seqID[-1]){
      points(y[, 1], y[, j+1], type = "l", col = ourcol[j%%8+1], lwd = lwd) 
    }
    abline(h = 0)
    abline(v = xi, lty = 3, lwd = 0.5)
  }
  # par(mfrow = c(1,1))

  layout(matrix(1:1, 1, 1))
}

#---------------------------------------------#
#--- Sec 11 AuxFun for handling 'der' SLOT ---#
#---------------------------------------------#

#Function for extracting the full der matrix given supp 
full_der = function(S, supp, n){
  # S: matrix, partial der matrix
  # supp: matrix, corresponding support
  # n: inner knots
  d = dim(supp)[1]
  rSS = matrix(0, nrow = n+2, ncol = dim(S)[2])
  SS = der_split(S, supp)
  
  for(i in 1:d){
    rSS[supp[i,1]:supp[i,2], ] = SS[[i]]
  }
  
  return(rSS)
}

# S der matrix
# function for extracting a sub der matrix corresponding to derivative as a list of matrices
# so that the positions of the derivatives for the components do not need to be calculated
der_split = function(S, supp){
  d = dim(supp)[1]
  dd = supp[,2]-supp[,1] 
  st = 1 # start row of der
  SS = list()
  for(i in 1:d){
    en = st + dd[i]
    SS[[i]] = S[st:en, ]
    st = en+1
  }
  return(SS)
}

#The following transforms the symmetrically around center matrix representation of the derivatives 
#to the one-sided representation (RHS-limits for the highest derivative)
sym2one_single = function(S, inv=FALSE){
  n = dim(S)[1]
  l = dim(S)[2]
  d = ceiling(n/2)
  if(inv){
    S[d:n, l] = c(0, S[d:(n-1), l])
    if(n/2-floor(n/2) == 0){S[d,l] = S[d+1,l]}
  }else{
    d = d+1
    S[(d-1):n, l] = c(S[d:n, l], 0)
  }
  return(S)
}

#----------------------------------------------#
#--- Sec 12 AuxFun for handling 'supp' SLOT ---#
#----------------------------------------------#

# function for checking if the support is a full range support
# supp: support matrix
# return logic, TRUE->full range suppost;
is.full = function(supp, full_supp){
  res = supp[1,1] == full_supp[1,1] & supp[1,2] == full_supp[1,2]
  return(res)
}

# function for taking union of two supports
# assist for function 'supp_union' and 'linear_engine_reduced'
# supp1: support matrix d1 X 2
# supp2: support matrix d2 X 2
# return: the union of supp1 and supp2
supp_union_bin = function(supp1, supp2){
  # internal knots within support
  id_1 = NULL
  for(i in 1:dim(supp1)[1]){
    id_1=c(id_1,supp1[i,1]:supp1[i,2])
  } #the indexes in the support
  id_1 = setdiff(id_1, c(supp1))
  id_2 = NULL
  for(i in 1:dim(supp2)[1]){
    id_2=c(id_2,supp2[i,1]:supp2[i,2])
  } #the indexes in the support
  id_2 = setdiff(id_2, c(supp2))
  # union of internal knots
  temp = sort(union(id_1, id_2))
  # separate temp if the differenc is larger than 2
  # record the number of knots within each parts
  num = diff(c(0, which(diff(temp) != 1), length(temp)))
  # generate the support
  supp = matrix(0, nrow = length(num), ncol = 2)
  for(i in 1:length(num)){
    temp_id = head(temp, num[i])
    supp[i,1] = head(temp_id,1)-1
    supp[i,2] = tail(temp_id,1)+1
    temp = temp[-(1:num[i])]
  }
  # combine two intervals shaping one knots
  temp = table(sort(c(supp)))
  supp = matrix(as.numeric(names(temp))[temp==1], byrow = TRUE, ncol = 2)
  return(supp)
}

#----------------------------------------------#
#------- Sec 13 AuxFun for 'project()' --------#
#----------------------------------------------#

#Computing inner product between data treated as a piecewise constant function and 
#the orthonormal spline basis
# @param fdata m x (N+1) matrix representing a sample of size N of descretized functional data; 
# In the first column, there are m values of the argument and in the remaining columns the values corresponding to
# these arguments. It is assumed that the argument are in the increasing order. 
# @param 'Splinet' object, made of n-k-1 splines that have one component of the support each; where n is the number
# of knots (including the endpoints), k is the order of smoothness. 
# It is intendent for 'Splinets'-object that represent bases but will work for any 'Splinet'-object as long as it 
# is made of splines with singelton component support.
# @return N x (n-k-1) matrix of the inner products of piecewise constant functions made of the fdata and 
# the splines in the splinet object. 

innerdb = function(fdata,basis){
  
  xi = basis@knots
  k = basis@smorder
  supp = basis@supp #it must be a 1 x 2 matrix form, not a vector
  S = basis@der
  n_basis=length(S)
  newS = S #the list to keep the new entries of the matrices of the derivatives for integral
  
  taylor = basis@taylor
  d = length(S)
  n = length(xi)-2 #the number of internal knots (without endpoints)
  
  N=dim(fdata)[2]-1
  x=fdata[,1,drop=FALSE]
  y=fdata[,-1,drop=FALSE]
  
  #the matrix for the result:
 # A=matrix(0,N,n-k+1) #n here counts only the internal knots
  A=matrix(0,N,n_basis) #n here counts only the internal knots
 
  
  #####
  #the main part
  ####
  #for(i in 1:(n-k+1)){#looping through the elements of the basis
  for(i in 1:(n_basis)){#looping through the elements of the basis 
    #the discretized data are treated as a function that between the two subsequentive arguments takes the RHS value
    rangeknots=xi[c(supp[[i]][1,1],supp[[i]][1,2])]
    ind_in_range=(x<rangeknots[2]) & (x>=rangeknots[1]) #indexes in the arguments that the values of the piecewise 
    xx=matrix(c(x[ind_in_range,],xi[supp[[i]][1,2]]),ncol=1)
    
    
    #constant are relevant for computation of the inner product with the given element of the basis
    loc_spl=subsample(basis,i) #extracting the basis element and restricting it to the support
    loc_spl@knots=loc_spl@knots[loc_spl@supp[[1]][1,1]:loc_spl@supp[[1]][1,2]] #restricting knots
    if(loc_spl@equid == FALSE){
    loc_spl@taylor=loc_spl@taylor[loc_spl@supp[[1]][1,1]:(loc_spl@supp[[1]][1,2]-1),] #adjusting Taylor
    }
    loc_spl@supp[[1]]=matrix(c(1,loc_spl@supp[[1]][1,2]-loc_spl@supp[[1]][1,1]+1),1,2) #adjusting support
    loc_spl@type='sp'
    #the spline is restricted to the support so the indefinite integral will be computed over relevant range
    
    intgr_os=integra(loc_spl)
    
    intgr_val=evspline(intgr_os,x=xx)
    intgr_val=intgr_val[,-1,drop=FALSE]
    diff_int=diff(intgr_val)
    
    A[,i]=t(y[ind_in_range,,drop=FALSE])%*%diff_int
  }
  
  
  return(A)
}

# Old one, working but less efficient
# #Computing inner product between data treated as a piecewise constant function and 
# #the orthonormal spline basis
# innerdb=function(fdata,basis){
#   
#   intgr_os=integra(basis) #the splines being integrals of the elements of the basis (spline of the order increased by one that 
#   #do not satisfy the boundary condition on the RHS end)
#   
#   A=evspline(intgr_os,x=c(fdata[,1],tail(basis@knots,n=1))) #evaluating the inner products in the embedding space
#   A=A[,-1]
#   DA=diff(A) #This matrix is new_n x n_basis, where new_n is the number of rows in new_fdsp
#   #n_basis=n_knots - k -1 is the number of elements in the basis
#   
#   A=t(fdata[,-1,drop=FALSE])%*%DA #This is N x n_basis, it represents the inner products of 
#   #the input treated as the piecewise constant function with the 
#   #elements of the basis. 
#   
#   return(A)
# }


#####

#----------------------------------------------#
#----- Sec 14 Auxfun for periodic Splinet -----#
#----------------------------------------------#



# function used in periodic splinets for generating a small interval cosists of the last k knots
# of the set of knots xi and k extra knots that preserves the distanse of the first k knots of 
# the set xi.
# @param xi numeric \code{n+2} vector, a vector of n+2 knot locations presented in the increasing order and without ties.
# @param k non-negative integer, the smoothnes order of the splines, i.e. the highest order of non-zero derivative. 

# @return numeric \code{2k} vector, a vector of 2k knot consists of the first and the last k knots of xi.

# @return numeric \code{2k} vector, a vector of 2k knot consists of the first and the last k knots of xi.

extraknots=function(xi, k = 3){
  diff=diff(xi)
  l=length(xi)
  xi=xi[l]
  for (i in 1:k) {
    xi=append(xi, max(xi)+ diff[i])
    xi=append(min(xi)-diff[l-i],xi)
  }
  return(xi)
  
}


# plot splines in polar coordinate
# 
# First, two polar coordingat plotting functions mimicking those 
#in the cartesian coordinates
polar_plot=function(theta,r,asp = 1,axes=F,xlab='',ylab='',type='l',zoom=1, ...){
  x=r*cos(theta)
  y=r*sin(theta)
  M=max(r)
  plot(x,y,asp = 1,axes=F,xlab='',ylab='',type=type, xlim=c(-M,M),ylim=c(-zoom*M,zoom*M), ...)
}
polar_lines=function(theta,r, ...){
  x=r*cos(theta)
  y=r*sin(theta)
  lines(x,y, ...)
}
#

#######################
# x - represents the points on interval [0,1] where splines will be evaluated for the plots
# sID - the indicies of the input splines to be plotted

plot.spline.p=function(object, x = NULL, sID = NULL, vknots=T, ...){

  if(length(sID) == 0){#if the indices of the splines are not specified
    sID=1:length(object@der)
  }
  Nsp=length(sID)
  n_so = length(object@der)
  xi = object@knots
  k = object@smorder
  n = length(xi)-2
  if(length(x)==0){
    x= (1/360)*seq(0,360, by = 1) #Change from degrees to interval [0,1]
  }
  
  y = evspline(object,x, sID = sID)
  y=y[,-1,drop=F]
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  
  l=dim(y)[2]
  M=max(y)
  a=log(2)/M
  c=rep(2, length(t))
  
  
  #entering the first functional datum (we change the argument to radians)
  polar_plot( 2*pi*x, exp(a*y[,1]), type="l",  col = ourcol[1%%8+1], zoom=exp(a*M)/exp(a*max(y[,1])), ...)
  #also scaling to match the scale with the largest value in all data in y
  
  #entering the remaining functional data
  if(l>1){
    for (i in 2:l) {
      polar_lines(2*pi*x, exp(a*y[,i]), type="l",  col = ourcol[i%%8+1], ...)
    } 
  }
  
  ####### Auxilary plotting
  #Plotting circles at the zero level, the maximal level, and the minimal level
  polar_lines(c(0,0),c(0,2),lty=1,col='black')
  one=rep(1, length(x))
  two=rep(2, length(x))
  
  polar_lines(2*pi*x, one, lty=1,col="black") #The circle at the level zero
  graphics::text(1+0.1,0.1,labels=0)
  polar_lines(2*pi*x, two, lty=3,col="black") #The circle at the level M
  graphics::text(2.1,0.1,labels=signif(M,3))
  m=min(y)
  polar_lines(2*pi*x, rep(exp(a*m),length(x)), lty=3,col="black") #The circle at the level m
  graphics::text(exp(a*m)-0.1,0.1,labels=signif(m,2))
  
  #Plotting the dashed lines at the knot locations
  if(vknots==TRUE)
  {
    tt=2*pi*object@knots
    for(i in 1:length(tt)){
      polar_lines(c(tt[i],tt[i]),c(0,2),lty=3,col='black')
    } 
  }
  #the end of knot plotting
  ###### End Auxilary plotting
}


#######
## Plotting bases of periodic splines
#######################
# x - represents the points on interval [0,1] where splines will be evaluated for the plots
# sID - the indicies of the input splines to be plotted
######
plot.splinet.p = function(object, x = NULL, sID = NULL, vknots=TRUE, type='l', bty="n",col='deepskyblue4',
                          lty=1, lwd=2, xlim=NULL, ylim = NULL, xlab="", ylab = "", ...){
  
  if(length(x) == 1){
    stop("x should not be a single value") # x - represent values at which the splines are 
      #evaluated to create graphs  if NULL the equally spaced grid of 360 points is chosen.
  }
  if(length(x)==0){ #x==NULL case
    x= (1/360)*seq(0,360, by = 1) #Change from degrees to interval [0,1]
  }
  
  if(length(sID) == 0){#If the indices are not specified all elements are plotted 
    sID=1:length(object@der) #the number of all input spline
  }
  Nsp=length(sID)
  n_so = length(object@der)
  xi = object@knots
  k = object@smorder
  n = length(xi)-2
  y = evspline(object,x=x) #Here the values for plots are evaluated
  
  
  # plot setups
  LWD = 2
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  
  net_str = net_structure(n-k+1, k) #extracts the levels and tuplets and individual splines for the number of splines n-k+1 of order k
  n_level = max(net_str[, 2]) #the first column of 'net_str' is 1:(n-k+1), the second identifies the level index 
                              #(top is the lowest level, bottom the highest, which opposite which we normally denoted in the paper), 
                              #For periodic splines there is one more k-tuple at the end, which is not 
                              #considered here. 
  
  
  ma=max(y[,-1,drop=F]) #the largest value within the splines
  
  a=log(2)/ma #
  c=rep(2*(n_level), length(x)) #all vallues will be within this circle
  
  #The highest level first
  
  
  polar_plot(2*pi*x, c, lty= 3, ...) #maximal level circle for the outmost level 
                                  #(the with initial splines in the basis)
  for(i in 1:(n_level-1)){#maximal circles for the remaining levels (they will be also -infinity levels for one up levels)
    c=rep(2*(n_level-i), length(x))
    polar_lines(2*pi*x, c, lty= 3) 
  }
  
  #The circles corresponding to zero level at each pyramid level
  for(i in 1:(n_level)){#maximal circles for the remaining levels (they will be also -infinity levels for one up levels)
    c=rep(2*(n_level-i)+1, length(x))
    polar_lines(2*pi*x, c, lty= 3) 
  }
  
  #Plotting the dashed lines at the knot locations
  if(vknots==TRUE)
  {
    tt=2*pi*xi
    for(i in 1:length(tt)){
      polar_lines(c(tt[i],tt[i]),c(0,2*(n_level)),lty=3)
    } 
  }
  
  #Identification of the splines to be plot on each level
  #Level one: 
  for(l in 1:n_level){#Plotting on subsequent levels
    seqID = sID[which(net_str[,2] == n_level-l+1)]
  
    if(length(seqID)>1){
      for(j in seqID){
        polar_lines(2*pi*x,(2*(n_level-l))+ exp(a* y[,(j+1)]), type = "l",  col = ourcol[j%%8+1])
      }
    }
  }
  
  for(j in (n_so-k+2):(n_so+1)){
    polar_lines(2*pi*x,exp(a* y[,j]), type = "l",  col = ourcol[j%%8+1])
  }
  # par(mfrow = c(1,1))
  layout(matrix(1:1, 1, 1))
}

#####
#----------------------------------------------------#
#---- Sec 15 two main functions for fun_splinets ----#
#----------------------------------------------------#

# Function for generating B-splines and their orthogonalization
#' @keywords internal

# splinet1 was fun_splinet in Version 1 of Splinet package.

splinet1 = function(knots=NULL, smorder = 3, type = 'spnt', Bsplines=NULL, norm=F){
  
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
    #'Splinets'-class defintion
    so@type = "bs"
    supp = cbind(1:(n-k+1),(1:(n-k+1))+k+1)
    supp = so@supp = lapply(seq_len(nrow(supp)), function(i) supp[i,,drop=F]) #seq_len(n)=1:n, 
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
        SS[[i]] = sym2one(rbind(S, numeric(k+1)), supp[[i]], inv = T)
      }
    } else{
      for(i in 1:(n-k+1)){
        SS[[i]] = sym2one(rbind(S[(i+(i-1)*k):(i+i*k),], numeric(k+1)), supp[[i]], inv = T)
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
  a = array(1/sqrt(gramian(Bsplines, norm_only = T)), dim=c((n-k+1),1))
  so = Bsplines
  for(i in 1:length(a)){
    so@der[[i]] = a[i]*so@der[[i]]
  }
  if(norm==T){Bsplines=so} #Normalization of the output B-splines
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

# Function for generating periodic B-splines and their orthogonalization

#splinet2 function uses splinets1 function which was fun_splinet in version 1 of the package.

splinet2 = function(knots, smorder = 3, type = 'spnt', Bsplines=NULL, norm=F){
  #In the case knots are not in the increasing order they are sorted
  if(!is.null(knots)){
  if(min(diff(knots))<0){
    knots=sort(unique(knots))
    cat("Knots were not given in the strictly increasing order, which is required.\n
          Ordered  knots with removed ties are replacing the input values.\n")
  }
  }
  #
  
  if(!is.null(Bsplines)){ #inheriting the arguments if B-splines are in the input
   # stopifnot(is(Bsplines,"Splinets"))
    if(Bsplines@periodic!= TRUE){# no Splinets in the input thus knots have to be taken from the argument
       stop("The B-splines in the input are not periodic, which is required in the periodic case.")}
    else{
    knots=Bsplines@knots
    smorder=k=Bsplines@smorder
    n = length(knots) - 2
    Bsplines=Bsplines
    Bspl=Bsplines
    n_sop = length(Bsplines@supp)
    lBspl=length(Bspl@der)
    Bspl@supp=Bspl@supp[1:(lBspl-k)]
    Bspl@der=Bspl@der[1:(lBspl-k)]
    Bspl@periodic=FALSE
    }
  }else{
    #------------------------------#
    # S1: generating bspline basis #
    #------------------------------#
    k = smorder
    n = length(knots) -2
    sop = new("Splinets", knots = knots, smorder = k)
    sop@type = "bs"
    sop@periodic=TRUE
    
    #generating Bsplines over the given knots
    Bspl=splinet1(knots, k, type = 'bs') 
    Bspl=Bspl$bs
    S=Bspl@der
    supp=Bspl@supp
    lBspl=length(S)
    lM= dim(S[[1]])[1]
    
    #generating a small interval with 2k+1 knots
    eknots=extraknots(knots,k )
    
    #generating the k extra splines to complete the set of periodic splines
    espline=splinet1(eknots,k, type = 'bs' )
    eS=espline$bs@der
    esupp=espline$bs@supp
    lespline=length(eS)
    
    #transferring the der matrices from symmetric to one sided
    for (i in 1:k) {
      eS[[i]]= sym2one(eS[[i]], esupp[[i]])
      
    }
    
    #changing support of the extra splines 
    supp1 = cbind((lBspl+1):(lBspl+k),(n+2):(n+2))
    supp1 = lapply(seq_len(lespline), function(i) supp1[i,,drop=F])
    supp2 = cbind(1:1,2:(k+1))
    supp2 = lapply(seq_len(lespline), function(i) supp2[i,,drop=F])
    esupp = lapply(seq_len(lespline), function(i) rbind(supp2[[i]],supp1[[i]]))
    
    
    ####changing der matrix of the extra splines
    
    Mder1=lapply(seq_len(k), function(i) eS[[(i)]][(1):(lM-i),,drop=F])
    Mder2=lapply(seq_len(k), function(i) eS[[(i)]][(lM-i):(lM),,drop=F])
    
    # In this step we changed the last row in Mder1 to zeros and transferring
    #the der matrices from one sided to symmetric
    for (i in 1:k) { 
      D=dim(Mder1[[i]])
      Mder1[[i]][D[1],]=0 
      Mder1[[i]]=sym2one(Mder1[[i]],supp1[[i]],inv=T)
    }
    
    for (i in 1:k) { 
      D=dim(Mder2[[i]])
      Mder2[[i]]=sym2one(Mder2[[i]],supp2[[i]],inv=T)
    }
    
    
    Mder= lapply(seq_len(lespline), function(i) rbind(Mder2[[i]],Mder1[[i]]))
    
    for (i in 1:k) { ## extra k splines that have the end points (the first and the last points of the iterval connected)
      ##in its support
      S[[lBspl+i]]=Mder[[i]]
      supp[[lBspl+i]]=esupp[[i]]
    }
    
    
    
    sop@der = S
    sop@supp=supp
    n_sop = length(S)
    Bsplines=sop
  } #The end of generation of the periodic B-splines if they are not in the input
  #------------------------------#
  # S2: normalization of bspline basis #
  #------------------------------#  
  # normalization of the B-splines (one does not assume that the input splines are orthogonalized)
  a = array(1/sqrt(gramian(Bsplines, norm_only = T)), dim=c((n+1),1))
  so = Bsplines
  for(i in 1:length(a)){
    so@der[[i]] = a[i]*so@der[[i]]
  }
  if(norm==T){Bsplines=so} #Normalization of the output B-splines
  splnt=list("bs"=Bsplines) #The B-spline part of the output list
  n_so=length(so@der) #The number of osplines in the basis
  #------------------------------#
  # S3: Orthogonalization of B-splines #
  #------------------------------#  
  # S3.1: gram schmidt obasis
  if(type == 'gsob'){
    H = bandmatrix(knots, k, so@der, so@supp)
    P = grscho(diag(n_so), H) #the generic algorithm for GS orthogonalization
    so = lincomb(so, t(P))
    so@type = 'gsob'
    splnt$os=so
  }
  
  # S3.2: twosided obasis
  if(type == 'twob'){
    two= splinet1(knots, k, Bsplines=Bspl,norm = T, type = 'twob')
    n_two=length(two$bs@der)
    ctwo= so
    for (j in 1:n_two) {
      ctwo@der[[j]]= two$os@der[[j]]
      ctwo@supp[[j]]= two$os@supp[[j]]
    }
    gram=gramian(ctwo, norm_only = F)
    # matrix with negative inner products
    A= -1*(gram-diag(diag(gram),nrow=(n_two+k))) + diag(diag(gram),nrow=(n_two+k))
    lA=dim(A)[1]
    A=A[(lA-k+1):lA,]
    
    d=dim(A)[2]
    
    for (i in 1:k) {
      for (j in (d-k+1):d) {
        if (!(j-i)==(d-k)) { A[i,j]=0}
      }
    }
    
    s = lincomb(ctwo,A)
    H=gramian(s)
    P = dyadicsplinet(A = diag(k), H = H, k = k)
    # generate the splinet os
    ss = lincomb(s, t(P))
    
    if (length(ss@supp)>0) {
      for (j in 1:k) {
       ctwo@der[[n_two+j]]= ss@der[[j]]
       ctwo@supp[[n_two+j]]= ss@supp[[j]]
      }
    }else{
      for (j in 1:k) {
        ctwo@der[[n_two+j]]= ss@der[[j]]
        ctwo@supp[[n_two+j]]= matrix(c(1,(n+2)),nrow=1)
      }
    }
    ctwo@type = 'twob'
    splnt$os=ctwo
  }
  
  
  # S2.3: splinet
  if(type == 'spnt'){
    #generating the original splinets  over the given knots
    net_splinet=splinet1(knots, k, Bsplines=Bspl,norm = T, type = 'spnt')
    n_os=length(net_splinet$os@der)
    
   
    # # choosing the first and the last tuplet at each level.
    M=net_splinet$os
    net_str = net_structure(n-k+1, k)
    n_level = max(net_str[, 2])
    
    
    M_net=list()
    M_supp=list()
    s=1
    seqID = net_str[which(net_str[,2] == 1),1 ]
    for (j in seqID) {
      M_net[[s]]= net_splinet$os@der[[j]]
      M_supp[[s]]= net_splinet$os@supp[[j]]
      s=s+1
    }
    
    for(i in 2:n_level) {
      mmax=max(net_str[which(net_str[,2] == i),3])
      mmin=min(net_str[which(net_str[,2] == i),3])
      for (j in c(mmin,mmax))
      {
        seqID=net_str[net_str[,2] == i & net_str[,3]==j, 1]
        for (l in seqID)
        {
          M_net[[s]]= net_splinet$os@der[[l]]
          M_supp[[s]]= net_splinet$os@supp[[l]]
          s=s+1
        }
      }
    }
    
    # adding the extra splines that need to be orthogonalized with respect to
    # the first and the last tuplet at each level.
    
    lM_net=length(M_net)
    for (i in 1:k) {
      M_net[[lM_net+i]]= so@der[[n_sop-k+i]]
      M_supp[[lM_net+i]]= so@supp[[n_sop-k+i]]
    }
    
    
    M@supp=M_supp
    M@der=M_net
    M@type='sp'
    gram=gramian(M, norm_only = F)
    
    
    # matrix with negative inner products
    A= -1*(gram-diag(diag(gram),nrow=(lM_net+k))) + diag(diag(gram),nrow=(lM_net+k))
    lA=dim(A)[1]
    A=A[(lA-k+1):lA,]
    
    d=dim(A)[2]
    
    for (i in 1:k) {
      for (j in (d-k+1):d) {
        if (!(j-i)==(d-k)) { A[i,j]=0}
      }
    }
    
    s = lincomb(M,A)
    H=gramian(s)
    P = dyadicsplinet(A = diag(k), H = H, k = k)
    # generate the splinet os
    ss = lincomb(s, t(P))
    
    
    
    
    
    if (length(ss@supp)>0) {
      for (j in 1:k) {
        net_splinet$os@der[[n_os+j]]= ss@der[[j]]
        net_splinet$os@supp[[n_os+j]]= ss@supp[[j]]
      }
    }
      else
      {
        for (j in 1:k) {
          net_splinet$os@der[[n_os+j]]= ss@der[[j]]
          net_splinet$os@supp[[n_os+j]]= matrix(c(1,(n+2)),nrow=1)
        }
      }
  
    splnt$os=net_splinet$os
    splnt$os@periodic=TRUE
    res=aug_bound(n_so,k)
    
    
    if(res$n_U == 0 & res$n_D == 0){
      splnt$os@type = "dspn"
    }else{
      splnt$os@type = "spnt"
    }
  }
  return(splnt)
}

