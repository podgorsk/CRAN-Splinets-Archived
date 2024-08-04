#-----------------------------------------------------#
#---------------------Subsampling---------------------#
#-----------------------------------------------------#

#Example with different support ranges, the 3rd order
n=25; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1; k=3
supp=list(t(c(2,12)),t(c(4,20)),t(c(6,25))) #defining support ranges for three splines

#Initial random matrices of the derivative for each spline
set.seed(5)
SS1=matrix(rnorm((supp[[1]][1,2]-supp[[1]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[[2]][1,2]-supp[[2]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[[3]][1,2]-supp[[3]][1,1]+1)*(k+1)),ncol=(k+1)) 

spl=construct(xi,k,SS1,supp[[1]]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[[2]],'CRFC')

#See 'gather' function for more details on what follows
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[[3]],'CRLC')
spl=gather(spl,nspl) #the third is added

#Replicating by subsampling with replacement
sz=length(spl@der)
ss=sample(1:sz,size=10,rep=TRUE)

spl=subsample(spl,ss)
is.splinets(spl)[[1]]

spl@supp
spl@der

#Subsampling without replacements
ss=c(3,8,1)
sspl=subsample(spl,ss)

sspl@supp
sspl@der

is.splinets(sspl)[[1]]

#A single spline sampled from a 'Splinets' object
is.splinets(subsample(sspl,1))

