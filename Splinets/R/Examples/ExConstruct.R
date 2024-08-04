#-------------------------------------------------------------#
#---Building 'Splinets' using different derviative matching---#
#-------------------------------------------------------------#
n=17; k=4
set.seed(5)
xi=sort(runif(n+2)); xi[1]=0; xi[n+1]=1 

#Random matrix of derivatives -- the noise (wild) case to be corrected
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) #construction of an object, the order of knots is corrected
is.splinets(spl)[[1]] #validation

spl=construct(xi,k,S,mthd='CRFC') #another method of the derivative matching
is.splinets(spl)[[1]]

spl=construct(xi,k,S,mthd='CRLC') #one more method
is.splinets(spl)[[1]]             

#-----------------------------------------------------#
#---------Building not over the full support----------#
#-----------------------------------------------------#
set.seed(5)
n=20; xi=sort(runif(n+2));xi[1]=0;xi[n+2]=1

spl=construct(xi,k,S) #construction of a spline as the 'Splinets' object over the entire range
is.splinets(spl)[[1]] #verification of the conditions

supp=c(3,17) #definition of the single interval support
SS=matrix(rnorm((supp[2]-supp[1]+1)*(k+1)),ncol=(k+1)) #The matrix of derivatives 
                                                       #over the support range
sspl=construct(xi,k,SS,supp=supp) #construction of a spline as the 'Splinets' object 
                                  #with the given support range
is.splinets(sspl)[[1]] #Verification 
sspl@knots
sspl@supp
sspl@der
