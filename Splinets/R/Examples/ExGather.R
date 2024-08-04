#-------------------------------------------------------------#
#-----------------Grouping into a 'Splinets' object-----------#
#-------------------------------------------------------------#
#Gathering three 'Splinets' objects using three different
#method to correct the derivative matrix
set.seed(5)
n=13;xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1; k=4
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) #constructing the first correct spline
spl=gather(spl,construct(xi,k,S,mthd='CRFC')) #the second and the first ones
spl=gather(spl,construct(xi,k,S,mthd='CRLC')) #the third is added

is.splinets(spl)[[1]] #diagnostic
spl@supp #the entire range for the support

#Example with different support ranges, the 3rd order
set.seed(5)
n=25; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1; k=3
supp=list(t(c(2,12)),t(c(4,20)),t(c(6,25))) #support ranges for three splines

#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[[1]][1,2]-supp[[1]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[[2]][1,2]-supp[[2]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[[3]][1,2]-supp[[3]][1,1]+1)*(k+1)),ncol=(k+1)) 

spl=construct(xi,k,SS1,supp[[1]]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[[2]])
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[[3]])
spl=gather(spl,nspl) #the third is added
is.splinets(spl)[[1]]

spl@supp
spl@der