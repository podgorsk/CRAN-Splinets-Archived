#-------------------------------------------------------#
#--- Generating the deriviative functions of splines ---#
#-------------------------------------------------------#
n=13; k=4
set.seed(5)
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
spl=construct(xi,k,matrix(rnorm((n+2)*(k+1)),ncol=(k+1))) #constructing three splines
spl=gather(spl, construct(xi,k,matrix(rnorm((n+2)*(k+1)),ncol=(k+1)))) 
spl=gather(spl, construct(xi,k,matrix(rnorm((n+2)*(k+1)),ncol=(k+1)))) 
# calculate the derivative of splines
dspl = deriva(spl)
plot(spl)
plot(dspl)

#----------------------------------------------#
#--- Examples with different support ranges ---#
#----------------------------------------------#

n=25; k=3
xi=seq(0,1,by=1/(n+1));
set.seed(5)
#Defining support ranges for three splines
supp=matrix(c(2,12,4,20,6,25),byrow=TRUE,ncol=2)
#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[1,2]-supp[1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[2,2]-supp[2,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[3,2]-supp[3,1]+1)*(k+1)),ncol=(k+1)) 
spl=construct(xi,k,SS1,supp[1,]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[2,])
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[3,])
spl=gather(spl,nspl) #the third is added

der_spl = deriva(spl)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
plot(der_spl)
plot(spl)
par(mfrow=c(1,1))
