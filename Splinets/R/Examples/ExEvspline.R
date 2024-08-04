#---------------------------------------------#
#-- Example piecewise polynomial vs. spline --#
#---------------------------------------------#
n=20; k=3; xi=sort(runif(n+2))
sp=new("Splinets",knots=xi) 

#Randomly assigning the derivatives -- a very 'wild' function.
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
sp@supp=list(t(c(1,n+2))); sp@smorder=k; sp@der[[1]]=S

y = evspline(sp)
plot(y,type = 'l',col='red')

#A correct spline object
nsp=is.splinets(sp)
sp2=nsp$robject 
y = evspline(sp2)
lines(y,type='l')

#---------------------------------------------#
#-- Example piecewise polynomial vs. spline --#
#---------------------------------------------#
#Gathering three 'Splinets' objects using three different
#method to correct the derivative matrix
n=17; k=4; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1

S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1)) # generate a random matrix S

spl=construct(xi,k,S) #constructing the first correct spline
spl=gather(spl,construct(xi,k,S,mthd='CRFC')) #the second and the first ones
spl=gather(spl,construct(xi,k,S,mthd='CRLC')) #the third is added
y = evspline(spl, sID= 1)
plot(y,type = 'l',col='red')

y = evspline(spl, sID = c(1,3))
plot(y[,1:2],type = 'l',col='red')
points(y[,c(1,3)],type = 'l',col='blue')

#sID = NULL
y = evspline(spl)
plot(y[,1:2],type = 'l',col='red',ylim=range(y[,2:4]))
points(y[,c(1,3)],type = 'l',col='blue')
points(y[,c(1,4)],type = 'l',col='green')

#---------------------------------------------#
#--- Example with different support ranges ---#
#---------------------------------------------#
n=25; k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
#Defining support ranges for three splines
supp=matrix(c(2,12,4,20,6,25),byrow=TRUE,ncol=2)
#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[1,2]-supp[1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[2,2]-supp[2,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[3,2]-supp[3,1]+1)*(k+1)),ncol=(k+1)) 
spl=construct(xi,k,SS1,supp[1,]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[2,],'CRFC')
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[3,],'CRLC')
spl=gather(spl,nspl) #the third is added

y = evspline(spl, sID= 1)
plot(y,type = 'l',col='red')

y = evspline(spl, sID = c(1,3))
plot(y[,1:2],type = 'l',col='red')
points(y[,c(1,3)],type = 'l',col='blue')

#sID = NULL -- all splines evaluated
y = evspline(spl)
plot(y[,c(1,3)],type = 'l',col='red',ylim=c(-1,1))
points(y[,1:2],type = 'l',col='blue')
points(y[,c(1,4)],type = 'l',col='green')

