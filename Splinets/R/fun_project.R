#' @title Projecting into spline spaces
#' 
#' @description  The projection of splines or functional data into the linear spline space spanned over a given set of knots.  
#' @param fdsp \code{Splinets}-object or a \code{n x (N+1)} matrix, a representation of \code{N} functions to be projected to 
#' the space spanned by a \code{Splinets}-basis over a specific set of knots;
#' If the parameter is a \code{Splinets}-object containing \code{N} splines, 
#' then it is orthogonally projected or represented in the basis that is specified by other parameters. 
#' If the paramater is a matrix, 
#' then it is treated as \code{N} piecewise constant functions with the arguments in the first column and the corresponding values of the functions in the remaining \code{N} columns. 
#' @param knots vector, the knots of the projection space, together with \code{smorder} fully characterizes the projection space; This parameter is overridden by the SLOT \code{basis@knots} of the \code{basis} input if this one is not \code{NULL}.
#' @param smorder integer, the order of smoothness of the projection space; This parameter is overridden by  the SLOT \code{basis@smorder} of the \code{basis} input if this one is not \code{NULL}.
#' @param periodic logical, a flag to indicate if B-splines will be of periodic type or not; In the case of periodic splines, the arguments of the input and the knots need to be within [0,1] or, otherwise, an error occurs and a message advising the recentering and rescaling data is shown. 
#' @param basis \code{Splinets}-object, the basis used for the representation of the projection of the input \code{fdsp};
#' @param type string, the choice of the basis in the projection space used only if the \code{basis}-parameter is not given; The following choices are available 
#' \itemize{
#'   \item \code{'bs'} for the unorthogonalized B-splines,
#'   \item \code{'spnt'} for the orthogonal splinet (the default),
#'   \item \code{'gsob'} for the Gramm-Schmidt (one-sided) OB-splines,
#'   \item \code{'twob'} for the two-sided OB-splines.
#'  } 
#' The default is \code{'spnt'}. 
#' @param graph logical, indicator if the illustrative plots are to be produced:
#'  \itemize{
#'   \item the splinet used in the projection(s) on the dyadic grid,
#'   \item the coefficients of the projection(s) on the dyadic grid,
#'   \item the input function(s),
#'   \item the projection(s).
#'  } ;
#' @return The value of the function is a list made of four elements
#' \itemize{
#'   \item \code{project$input} -- \code{fdsp}, when the input is a \code{Splinets}-object or a matrix with the first column in an increasing order,
#'   otherwise it is the input numeric matrix after ordering according to the first column, 
#'   \item \code{project$coeff} -- \code{N x (n-k+1)} matrix of the coefficients of representation of the projection of the input in the splinet basis,
#'   \item \code{project$basis} -- the spline basis,
#'   \item \code{projedt$sp} -- the \code{Splinets}-object containing the projected splines.
#'  } 
#' @details The obtained coefficients \eqn{\mathbf A = (a_{ji})} with respect to the basis allow to evaluate the splines \eqn{S_j} in the projection according to 
#' \deqn{
#' S_j=\sum_{i=1}^{n-k-1} a_{ji} OB_{i}, \,\, j=1,\dots, N,
#'  }
#' where \eqn{n} is the number of the knots (including the endpoints), \eqn{k} is the spline smoothness order, 
#' \eqn{N} is the number of the projected functions and  \eqn{OB_i}'s consitute the considered basis. 
#' The coefficient for the splinet basis are always evaluated and thus, for example, 
#' \code{PFD=project(FD,knots); ProjDataSplines=lincomb(PFD$coeff,PFD$basis)}
#' creates a \code{Splinets}-object made of the projections of the input functional data in \code{FD}.
#' If the input parameter \code{basis} is given, then the function utilizes this basis and does not 
#' need to build it. However, if \code{basis} is the B-spline basis, then the B-spline orthogonalization is performed anyway, 
#' thus the computational gain is smaller than in the case when \code{basis} is an orthogonal basis. 
#' 
#' @export
#' @seealso \code{\link{refine}} for embeding a \code{Splinets}-object into the space of splines with an extended set of knots; 
#' \code{\link{lincomb}} for evaluation of a linear combination of splines; 
#' \code{\link{splinet}} for obtaining the spline bases given the set of knots and the smootheness order; 
#' @inheritSection Splinets-class References
#' @example R/Examples/ExProject.R
#' @importFrom grDevices dev.new
#' @importFrom graphics abline layout matplot par points
#' @importFrom utils capture.output head tail
#' 

project = function(fdsp, knots=NULL, smorder=3, periodic=FALSE, basis=NULL, type='spnt', graph=FALSE){


proj=list(input=NULL,coeff=NULL,basis=NULL,sp=NULL) #this is where output is kept

############################
###Checking Input values####
############################
 #In the first part, the proper input values are check so that 
 #the projection space is properly decided for.
  if(periodic==TRUE){
    #Checking if all arguments of functions are within [0,1]
    if(is.matrix(fdsp)){
      Mi=min(fdsp[,1])
      Ma=max(fdsp[,1])
      if(Mi<0||Ma>1){stop("The argument of the input data is not normalized to be in [0,1], which is required for the periodic case.\n
Consider renormalization.")
        }
    }else{
      Mi=min(fdsp$knots)
      Ma=max(fdsp$knots)
      if(Mi<0||Ma>1){stop("The knots of the input spline data are not normalized to be in [0,1], which is required for the periodic case.\n
Consider renormalization.")
        }
    }
    if(!is.null(knots)){
      Mi=min(knots)
      Ma=max(knots)
      if(Mi<0||Ma>1){stop("The specified knots are not normalized to be in [0,1], which is required for the periodic case.\n
Consider renormalization.")
      }
    }
      if(!is.null(basis)){
        Mi=min(basis$knots)
        Ma=max(basis$knots)
        if(Mi<0||Ma>1){stop("The knots of the input basis are not normalized to be in [0,1], which is required for the periodic case\n
Consider renormalization.")
        }
      }   
  }
  if(!is.null(knots)){
    if(min(diff(knots))<0){
      knots=sort(unique(knots))
      cat("The knots were not given in the strictly increasing order, which is required. The ordered knots with removed ties are replacing the input knot values.\n")
    }
  } #ordering knots and remove their ties (if present in the input) and transfering to the radian scale for the periodic case

  
  
  if(is.matrix(fdsp)){#in the case of the matrix input, we order the first column
  ld=length(fdsp[,1])
  isord=min(fdsp[2:ld,1]-fdsp[1:(ld-1),1])
  if(isord < 0){
    or=order(fdsp[,1])
    fdsp=fdsp[or,] #The ordering of the data with respect to the first column
    cat("The input numerical data were not given in the increasing order.\n The data are ordered and returned in the ordered form as the first value in the output list.\n")
  }
  
  proj$input=fdsp #The input data are returned as the first element in the output list.
  
  if(is.null(knots) & is.null(basis)){# no basis and knots in the input thus knots have to be taken from the argument
    knots=fdsp[,1]
    cat("The knots set to the first column in the input matrix.\n") #note that for periodic case these are now normalized.
  }#the first column of the input matrix is a natural choice for knots if they are not given 
}
  if(!is.null(basis)){#if the basis is provided then all parameters are taken from it
    
    if(basis@type=='sp'){stop("The 'basis' input has SLOT 'type='sp'' and thus is not a spline basis.")
      }else{
             knots=basis@knots 
             smorder=basis@smorder
             proj$basis=basis
            if(!is.matrix(fdsp)){#the input Splinets-object must agree in order with the basis
                  if(fdsp@smorder != smorder){stop("The order of the basis and the input do not agree.")}
             }
      }
    
  }else{#if the basis is not provided, then knots and smorder have to be given and the latter has to agree with the input
    if(!is.matrix(fdsp)){#the input Splinets-object overwrites smorder if they do not agree
      if(fdsp@smorder != smorder){
        cat("The order of the input which is ",fdsp@smorder," overwrites smorder=",smorder)
        smorder=fdsp@smorder
      }
      if(is.null(knots)){knots=fdsp@knots}#The knots are taken from the input if not given.
    }
  }
  
  #transformed to the radian-scale.
  k=smorder #for convenience
 #At this point all the inputs are resolved to be coherent
  
  ####################################
  ######The main part of the code#####
  ####################################
  
    if(is.matrix(fdsp)){ #the case of a discretized data input
      
      #The part of the functional data beyond the range of the projection knots will be neglected in 
      #the computaions below
      new_fdsp=fdsp
      ir=range(fdsp[,1])
      kr=range(knots)
      if((ir[1]<kr[1]) | ir[2]>kr[2]){
        drop=sum(fdsp[,1]<kr[1])+sum(fdsp[,1]>kr[2])
       cat("The range of the input data is larger than the range of knots in the projection space.\n
The ", drop, " values at the arguments outside the projection range will not affect the projection.\n")
        ind=(fdsp[,1]>=kr[1] & fdsp[,1]<kr[2]) #the indices of the input within the projection range
        new_fdsp=fdsp[ind,] #for the purpose of the computations extracting the data withing the projection knot range.
      }
     
      
      if(is.null(proj$basis)){#the basis needs to be built
        bsp=splinet(knots=knots, smorder = smorder, type = type, Bsplines=basis, periodic= periodic)
         #This in the orthogonal basis case builds both the B-splines and the OB-splines
        
        if(type =='bs'){# the case of the B-splines (non-orthogonal basis)
          
          proj$basis=bsp$bs #this basis is non-normalized
          
          #At the moment the basis is not normalized because the default 
          #in 'splinet' generates unnormalized version of the B-splines
          obs=splinet(Bsplines = bsp$bs,norm = TRUE, periodic = periodic) #computing the orthogonal splines and the input B-splines are normalized
          #the first argument is the  B-spline basis unnormalized         
          
          #Computing the coefficients of the projection
          N=diag(1/sqrt(gramian(bsp$bs,norm_only=TRUE))) #The normalizing matrix, It will be needed to express the coefficient 
                                                      #in the non-normalized B-splines
          
          H=bandmatrix(knots,k,obs$bs@der,obs$bs@supp) #bandmatrix() requires that the splines are normalized and the output
                                                       #obs from splinet() is normalized
          
          A=innerdb(new_fdsp,obs$os) #To be efficiently implemented 
          # #Computing the inner product of the piecewise constant functions with the orthonormal basis elements
          #   intgr_os=integra(obs$os) #the splines being integrals of the elements of the basis (spline of the order increased by one that 
          #                            #do not satisfy the boundary condition on the RHS end)
          #   
          #   A=evspline(intgr_os,x=c(new_fdsp[,1],tail(knots,n=1))) #evaluating the inner products in the embedding space
          #   A=A[,-1]
          #   DA=diff(A) #This matrix is new_n x n_basis, where new_n is the number of rows in new_fdsp
          #              #n_basis=n_knots - k -1 is the number of elements in the basis
          #   
          #   A=t(new_fdsp[,-1,drop=FALSE])%*%DA #This is N x n_basis, it represents the inner products of 
          #                                  #the input treated as the piecewise constant function with the 
          #                                  #elements of the basis. 
          # #The end of computing the inner products of the ON basis with the pieciwise constant functions.
          #   
          P=t(dyadiag(H,k,obs$bs@equid)) #the main band-matrix diagonalization algorithm to obtain the coefficients in the B-splines
          #it is designed so that the transpose of it represents the matrix of base change from the B-splinets
          #to the splinet
          proj$coeff=A%*%P%*%N #these are the coefficients in the unnormalized B-splines
        }else{ #orthonormal basis
          proj$basis=bsp$os
          A=innerdb(new_fdsp,proj$basis) #To be efficiently implemented 
          # #Computing the inner product of the piecewise constant functions with the orthonormal basis elements
          # intgr_os=integra(proj$basis) #the splines being integrals of the elements of the basis (spline of the order increased by one that 
          # #do not satisfy the boundary condition on the RHS end)
          # 
          # A=evspline(intgr_os,x=c(new_fdsp[,1],tail(knots,n=1))) #evaluating the inner products in the embedding space
          # A=A[,-1]
          # DA=diff(A) #This matrix is new_n x n_basis, where new_n is the number of rows in new_fdsp
          # #n_basis=n_knots - k -1 is the number of elements in the basis
          # 
          # A=t(new_fdsp[,-1,drop=FALSE])%*%DA #This is N x n_basis, it represents the inner products of 
          # #the input treated as the piecewise constant function with the 
          # #elements of the basis. 
          # #The end of computing the inner products of the ON basis with the piecewise constant functions.
          # 
            proj$coeff=A
        }
        #the end of the case when the basis has to be built 
      }else{#the basis is provided and assigned to 'proj$basis', 'knots' and 'smorder' have been already assigned
        
        type=proj$basis@type
        if(type=='bs'){
          # the case of the B-splines (non-orthogonal basis)
          #At the moment the basis is not normalized because the default 
          #in 'splinet' generates unnormalized version of the B-splines
          obs=splinet(Bsplines = proj$basis,norm = TRUE,periodic = periodic) #computing the orthogonal splines and the input B-splines is normalized
          
          #Computing the coefficients of the projection
          N=diag(1/sqrt(gramian(proj$basis,norm_only=TRUE))) #The normalizing matrix, can be the identity if th
          #the original basis is already normalized. It will be needed to express the coefficient 
          #in the non-normalized B-splines
          H=bandmatrix(knots,smorder,obs$bs@der,obs$bs@supp) #bandmatrix() requires that the splines are normalized
          
          A=innerdb(new_fdsp,obs$os) #To be efficiently implemented 
          
          # #Computing the inner product of the piecewise constant functions with the orthonormal basis elements
          # intgr_os=integra(obs$os) #the splines being integrals of the elements of the basis (spline of the order increased by one that 
          # #do not satisfy the boundary condition on the RHS end)
          # 
          # A=evspline(intgr_os,x=c(new_fdsp[,1],tail(knots,n=1))) #evaluating the inner products in the embedding space
          # A=A[,-1]
          # DA=diff(A) #This matrix is new_n x n_basis, where new_n is the number of rows in new_fdsp
          # #n_basis=n_knots - k -1 is the number of elements in the basis
          # 
          # A=t(new_fdsp[,-1,drop=FALSE])%*%DA #This is N x n_basis, it represents the inner products of 
          # #the input treated as the piecewise constant function with the 
          # #elements of the basis. 
          # #The end of computing the inner products of the ON basis with the pieciwise constant functions.
          
          P=t(dyadiag(H,k,obs$bs@equid)) #the main band-matrix diagonalization algorithm to obtain the coefficients in the B-splines
          #it is designed so that the transpose of it represents the matrix of base change from the B-splinets
          #to the splinet
          proj$coeff=A%*%P%*%N #these are the coefficients in the unnormalized B-splines
        }else{ #orthonormal basis
          
          #Computing the inner product of the piecewise constant functions with the orthonormal basis elements
          intgr_os=integra(proj$basis) #the splines being integrals of the elements of the basis (spline of the order increased by one that 
          #do not satisfy the boundary condition on the RHS end)
          
          A=evspline(intgr_os,x=c(new_fdsp[,1],tail(knots,n=1))) #evaluating the inner products in the embedding space
          A=A[,-1]
          DA=diff(A) #This matrix is new_n x n_basis, where new_n is the number of rows in new_fdsp
          #n_basis=n_knots - k -1 is the number of elements in the basis
          
          A=t(new_fdsp[,-1,drop=FALSE])%*%DA #This is N x n_basis, it represents the inner products of 
          #the input treated as the piecewise constant function with the 
          #elements of the basis. 
          #The end of computing the inner products of the ON basis with the piecewise constant functions.
          
          proj$coeff=A
          #the end of the orthonormal basis case
          }
        #the end of the case that the basis was given for the discrete input case    
        }
      #the end of the discrete input case  
    }else{#the case when the input is a splinet
      
      if(is.null(proj$basis)){#the basis needs to be built
        
        bsp=splinet(knots,smorder,type,periodic = periodic) #This in the orthogonal basis builds both the B-splines and the OB-splines

        if(type =='bs'){# the case of the B-splines (non-orthogonal basis)
          
          proj$basis=bsp$bs #the first argument is the  B-spline basis (non-normalized)          
          
          #At the moment the basis is not normalized because the default 
          #in 'splinet' generates unnormalized version of the B-splines
          obs=splinet(Bsplines = bsp$bs,norm = TRUE,periodic = periodic) #computing the orthogonal splines and the input B-splines are normalized
          
          #Computing the coefficients of the projection
          N=diag(1/sqrt(gramian(bsp$bs,norm_only=TRUE))) #The normalizing matrix will be needed to express the result in the n
          #non-normalized version
          
          H=bandmatrix(knots,smorder,obs$bs@der,obs$bs@supp) #bandmatrix() requires that the splines are normalized
          lnk=min(length(knots),length(fdsp@knots))
          unk=max(length(knots),length(fdsp@knots))
          
          if(lnk==unk & prod(knots[1:lnk]==fdsp@knots[1:lnk])){
          #the projection becomes just the representation case, the same sets of knots in the input and the projection space
            A=gramian(fdsp,Sp2=obs$os) #the coefficients of the expansion by the splinet
          }else{#the case of different sets of knots in the projection space vs. in the input
            embfdsp=refine(fdsp,newknots=knots) #embedding to a larger common spline space
            embbase=refine(obs$os,newknots=fdsp@knots)
            A=gramian(embfdsp,Sp2=embbase) #evaluating the inner products in the embedding space
          }
          P=t(dyadiag(H,k,obs$bs@equid)) #the main band-matrix diagonalization algorithm to obtain the coefficients in the B-splines
                                      #it is designed so that the transpose of it represents the matrix of base change from the B-splinets
                                      #to the splinet
          proj$coeff=A%*%P%*%N #these are the coefficients in the unnormalized B-splines
        }else{ #orthonormal basis
            proj$basis=bsp$os
            lnk=min(length(knots),length(fdsp@knots))
            unk=max(length(knots),length(fdsp@knots))
            
            if(lnk==unk & prod(knots[1:lnk]==fdsp@knots[1:lnk])){
              #the projection becomes just the representation case, the same sets of knots in the input and the projection space
               proj$coeff=gramian(Sp=fdsp,Sp2=bsp$os)
            }else{#the case of projecting to another space of splines
              embfdsp=refine(fdsp,newknots=knots) #embedding to a larger common spline space
              embbase=refine(bsp$os,newknots=fdsp@knots)
              proj$coeff=gramian(Sp=embfdsp,Sp2=embbase)  #evaluating the inner products in the embedding space
            }
            
          }
       #the end of the case when the basis is built 
      }else{#the basis is provided and assigned to 'proj$basis', 'knots' and 'smorder' have been already assigned
          
          type=proj$basis@type #The type is inherited from the type of the given basis
          if(type=='bs'){
            # the case of the B-splines (non-orthogonal basis)
            #At the moment the basis can be not normalized          
            
            #The splinet corresponding to the B-spline basis still needs to be computed
            obs=splinet(Bsplines = proj$basis,norm = TRUE,periodic = periodic) #The output will have B-splines orthogonalized
            
            #Computing the coefficients of the projection
            N=diag(1/sqrt(gramian(proj$basis,norm_only=TRUE))) #The normalizing matrix (could be identity if the original matrix
            #is already normalized)
            
            
            H=bandmatrix(knots,smorder,obs$bs@der,obs$bs@supp) #bandmatrix() requires that the splines are normalized
            lnk=min(length(knots),length(fdsp@knots))
            unk=max(length(knots),length(fdsp@knots))
            
            if(lnk==unk & prod(knots[1:lnk]==fdsp@knots[1:lnk])){
              #the projection becomes just the representation case, the same sets of knots in the input and the projection space
              A=gramian(fdsp,Sp2=obs$os) #the coefficients of the expansion by the splinet
            }else{#the case of different sets of knots in the projection space vs. in the input
              embfdsp=refine(fdsp,newknots=knots) #embedding to a larger common spline space
              embbase=refine(obs$os,newknots=fdsp@knots)
              A=gramian(embfdsp,Sp2=embbase) #evaluating the inner products in the embedding space
            }
            P=t(dyadiag(H,smorder,obs$bs@equid)) #the main band-matrix diagonalization algorithm to obtain the coefficients in the B-splines
            #it is designed so that the transpose of it represents the matrix of base change from the B-splinets
            #to the splinet
            proj$coeff=A%*%P%*%N #these are the coefficients in the unnormalized B-splines
          }else{ #orthonormal basis
            
            lnk=min(length(knots),length(fdsp@knots))
            unk=max(length(knots),length(fdsp@knots))
            
            if(lnk==unk & prod(knots[1:lnk]==fdsp@knots[1:lnk])){
              #the projection becomes just the representation case, the same sets of knots in the input and the projection space
              proj$coeff=gramian(Sp=fdsp,Sp2=proj$basis)
            }else{#the case of projecting to another space of splines
              embfdsp=refine(fdsp,newknots=knots) #embedding to a larger common spline space
              embbase=refine(proj$basis,newknots=fdsp@knots)
              proj$coeff=gramian(Sp=embfdsp,Sp2=embbase)  #evaluating the inner products in the embedding space
            }
          }
     #The end of the case that the basis was given for the splinet input case    
      }
      
    #The end of the splinet input case  
    }
  
proj$sp=lincomb(proj$basis,proj$coeff)

if(graph == TRUE){
  ourcol=c('deepskyblue4', 'darkorange3', 'goldenrod', 'darkorchid4',
           'darkolivegreen4', 'deepskyblue', 'red4', 'slateblue')
  
  
  #Plotting the basis - First graph
  dev.new()
  plot(proj$basis,main='The basis of the projection') #The plot method recognize periodic from non-periodic
  #The end of the first graph
  
  
  #Plotting the coefficents = Second graph
  dev.new() #A new graphical window
  if (!proj$basis@periodic) { #The regular (non-periodic) case
    n=length(proj$basis@knots) #The number of all knots (including the two endpoints)
    k=proj$basis@smorder #So that the number of basis elements is: n-k-1, here is the count: 
    #total number of knots is n, each B-spline extends over k+2 knots
    #if k is odd then we drop (k+1)/2 from each end side and the rest of knots
    #corresponds to the centers of the elements of the basis
    #if k is even we drop k/2 knots  from each end side and take the midpoints of the rest
    #as the centers of the elements of the basis (the number of knots . 
    if((k%%2)==0){
      x1=matrix(proj$basis@knots[(k/2+1):(n-k/2)],ncol=1)
      x=(x1[1:dim(x1)[1]-1,1]+x1[2:dim(x1)[1],1])/2 #these locations identify elements of the basis 
    }else{
      x=matrix(proj$basis@knots[((k+1)/2+1):(n-(k+1)/2)],ncol=1) #these locations identify elements of the basis
    }
  
  if(type=='spnt'){#special dyadic structure ploting for the splinet case
    n_so = length(proj$basis@der)
    xi = proj$basis@knots
    k = proj$basis@smorder
    n = length(xi) #the total number of knots (including endpoints)
    y = t(proj$coeff)
    xrange=range(xi)
    
    #Setting the plot parameters, margins, number of plots,
    mrgn=1.5
    
    par(mar = c(mrgn, 2*mrgn, mrgn, 2*mrgn))
    net_str = net_structure(n-k-1, k)
    n_level = max(net_str[, 2])
    layout(matrix(1:n_level, n_level, 1)) #rowwise n_level x 1 plots for each level
    
    #The top level is seperated to put the title at its top plot.
    seqID = net_str[which(net_str[,2] == 1), 1] #the second column in the net_str corresponds to the level.
    matplot(x[seqID],y[seqID,],xlim=xrange,pch='+',cex=2.5,col=ourcol,main='Coefficents of the  projections',xlab='',ylab='', bty="n")
    abline(h = 0)
    abline(v=x[seqID],lwd = 1.5,col = ourcol[seqID%%8+1])
    abline(v = xi, lty = 3, lwd = 0.5)
    if(n_level>1){
      for(i in 2:n_level){ #Going through the levels starting from the top and plotting the coefficients at the
        #center of corresponding element in the basis
        seqID = net_str[which(net_str[,2] == i), 1] #the second column in the net_str corresponds to the level.
        
        matplot(x[seqID],y[seqID,],xlim=xrange,pch='+',cex=2.5,col=ourcol,xlab='',ylab='', bty="n")
        abline(v=x[seqID], lwd = 1.5,col = ourcol[seqID%%8+1])
        abline(h = 0)
        abline(v = xi, lty = 3, lwd = 0.5)
      }
    }  
    
  # 
  }else{
    XY=cbind(x,t(proj$coeff))
    matplot(XY[,1],XY[,-1],pch='+',cex=1,col=ourcol,main='Coefficents of the projections',xlab='basis spline location',ylab='coefficient', bty="n")
  
    abline(h = 0)
  
    abline(v = x, lty = 3, lwd = 0.5)
  }
  }else{#the periodic case - ploting coefficients
    
    n=length(proj$basis@knots) #The number of all knots (including the endpoints so on the
    #circle first and the last are counted twice)
    k=proj$basis@smorder #So that the number of basis elements is: n-k-1 of the regular splines 
    #plus k that at added to cover the k basis elements corresponding knots around zero
    # These basis elements sit at the top level and come from the periodic boundary condition. 
    # The total number of the basis elements is the same as total number of knots,
    # i.e. n-1 (substraction of one comes from counting seperately the last and first knot).
    # See also non-periodic for the count.
    
    #The location on the circle where the values of the coefficient will be plot 
    if((k%%2)==0){
      x1=matrix(proj$basis@knots[(k/2+1):(n-k/2)],ncol=1) #The endpoints
      x=matrix((x1[1:dim(x1)[1]-1,1]+x1[2:dim(x1)[1],1])/2,ncol=1)  #The midpoints
    }else{
      x=matrix(proj$basis@knots[((k+1)/2+1):(n-(k+1)/2)],ncol=1)
    }
    #At this point x contains only the locations from the regular case, next, we add the knots
    #that cover the last k basis elements that come from the periodic boundary condition
    if((k%%2)==0){
      #Around the origin one needs to define the midpoints
      #First, k/2 midpoints below the origin 0=1, which requires the last k/2+1 knots
      #the last being treated as 1 not as 0
      x1=matrix(proj$basis@knots[(n-k/2):(n)],ncol=1)
      xx=matrix((x1[1:(k/2),]+x1[2:(k/2+1),])/2,ncol=1) #k/2 midpoints below 1=0
      x2=matrix(proj$basis@knots[1:((k/2)+1)],ncol=1) 
      xx=rbind(xx,matrix((x2[1:(k/2),]+x2[2:(k/2+1),])/2,ncol=1))
    }else{
      xx=matrix(c(proj$basis@knots[(n-(k-1)/2):(n-1)],proj$basis@knots[1:((k+1)/2)]),ncol=1) 
    }#starting from n-(k+1)/2+1 to (n-(k+1)/2+k)mod(n-1)  
    
    x=rbind(x,xx) #the locations are added at the end, at this point x is made of n-1 elements
   
      y = t(proj$coeff) #This is where the coefficients are kept
      
      #The midpoints in the polar coordinates. The scaling is in fact not needed because the konts are in [0,1]
      x= 2*pi*(x-proj$basis@knots[1])/(proj$basis@knots[n]-proj$basis@knots[1]) # transfer knots to radians
    
      ma=max(abs(y)) #range of the coefficients
      
      a=log(2)/ma #The scaling see also description of the plots for periodic splines and splinets
      
      #The knots endpoints in the polar coordiantes (the scaling is not needed because knots are suppose to be on [0,1])
      xxx=2*pi*(proj$basis@knots-min(proj$basis@knots))/(max(proj$basis@knots)-min(proj$basis@knots))
      
    if(type=='spnt'){#special dyadic structure ploting for the splinet case
        net_str = net_structure(n-k-1, k) #extracts the levels and tuplets and individual splines for 
        #the number of splines n-k-1 of order k (here n is the number of all knots including endpoints) 
        #the first column of 'net_str' are the numbers of knots 1:(n-k-1), the second identifies the level index 
        #(top is the lowest level, bottom the highest, which opposite which we normally denoted in the paper), 
        #this function is used from the zero boundary splines so the last k-tuple of the periodic case is not considered
        
        n_level = max(net_str[, 2])
        
        ##Coefficients of projection - setting coordinates for dyadic structure
        c=rep(2*n_level, length(xxx)) #all vallues will be within this circle
        polar_plot(xxx, c, lty= 3,  col = "red",  main = 'Coefficients of the projections')
        
        for(i in 1:(n_level-1)){#maximal circles for the remaining levels (they will be also -infinity levels for one up levels)
          c=rep(2*(n_level-i), length(xxx))
          polar_lines(xxx, c, lty= 3,  col = "red") 
          }
        #The circles corresponding to zero level at each pyramid level
        for(i in 1:(n_level)){#maximal circles for the remaining levels (they will be also -infinity levels for one up levels)
          c=rep(2*(n_level-i)+1, length(xxx))
          polar_lines(xxx, c, lty= 1) 
          }
        #Plotting the dashed lines at the knot locations
        
        for(i in 1:length(xxx)){
          polar_lines(c(xxx[i],xxx[i]),c(0,2*(n_level)),lty=3)
          } 

        
        for(i in 1:n_level){
          seqID = net_str[which(net_str[,2] ==  n_level-i+1), 1]
          for(j in seqID){ #the center of a spline and thus its position is shifted by k-1 forward from the position of the coefficient in y
            polar_lines(x[j],2*(n_level-i)+ exp(a* y[j]) ,type = "p",pch='+',cex=1.5,col=ourcol[j%%8+1])
            polar_lines(c(x[j],x[j]),c(2*(n_level-i)+1, 2*(n_level-i)+exp(a*y[j])), type = "l",lwd=2,col=ourcol[j%%8+1]) 
            }
          }
        for(j in (n-k):(n-1)){ #There is total n-1 coefficient to be plotted (this is the number of basis elements)
          polar_lines(x[j],exp(a* y[j]),type = "p",pch='+',cex=1.5,col=ourcol[(j+1)%%8+1])
          polar_lines(c(x[j],x[j]),c(1,exp(a* y[j])), type = "l",lwd=2,col=ourcol[(j+1)%%8+1])
          }
      }else{# Regular coefficents on a single circular plot
        c=rep(2, length(xxx)) #the maximal value of the radius in the polar represantion
        
        polar_plot(xxx, c, lty= 3,  col = "red",  main = 'Coefficients of the projections')
        c_1=rep(1, length(xxx))
        polar_lines(xxx, c_1, type = "l",col="blue")
        c_2=rep(exp(-a*ma), length(xxx)) #the lower bound of the radius in the polar represantion
        polar_lines(xxx, c_2, type = "l",lty=3, col="red")
        for(i in 1:length(x)){
          polar_lines(c(x[i],x[i]),c(0,2),lty=3,col='black')
        } 
        
        for (i in 1:length(y)) {
          polar_lines(x[i], exp(a*y[i]),type = "p",pch='+', cex=1.5, col = ourcol[i%%8+1])
          polar_lines(c(x[i],x[i]),c(1,exp(a*y[i])), type = "l",lwd=2,col=ourcol[i%%8+1]) 
        }
        graphics::text(exp(a*ma)-0.1,0.1,labels=signif(ma,2))
        graphics::text(1-0.1,0.1,labels=0)
        graphics::text(exp(-a*ma)-0.1,0.1,labels=signif(-ma,2))
      }
      
  }#The end of the second graphical illustration
  

 #The third graph - the input data
  dev.new()
  
  if(is.matrix(fdsp)){#The discrete data as the input
    if (!proj$basis@periodic) {# The non-periodic case
      matplot(fdsp[,1],fdsp[,-1],pch='.',cex=3,col=ourcol,main='Discrete functional data',xlab='',ylab='', bty="n")
      }else{
        #The periodic case
        x= 2*pi*fdsp[,1]
        l=dim(fdsp[,-1,drop=F])[2]
        ma=max(fdsp[,-1,drop=F])
        mi=min(fdsp[,-1,drop=F])
        a=log(2)/ma
        
        xxx=2*pi*seq(0,1,by=0.01)
        c=rep(2, length(xxx))
        polar_plot(xxx, c, lty= 3,  col = "red",  main = 'Discretized input data')
        c_1=rep(1, length(xxx))
        polar_lines(xxx, c_1, type = "l",col="blue")
        c_2=rep(exp(-a*mi), length(xxx)) #the lower bound of the radius in the polar represantion
        polar_lines(xxx, c_2, type = "l",lty=3, col="red")
        polar_lines(x, exp(a*fdsp[,-1,drop=F]),type = "p",pch='.', cex=3.5) #ploting data
        graphics::text(exp(a*ma)+0.1,0.1,labels=signif(ma,2))
        graphics::text(1+0.1,0.1,labels=0)
        graphics::text(exp(a*mi)-0.1,0.1,labels=signif(mi,2))
      }#End of the discrete data case
   }else{#The functional data as the input
    plot(fdsp,main='Input: functional data')
   }
  dev.new()
  #Plotting the projections
  plot(proj$sp, main='Data projection')
}

return(proj)
}