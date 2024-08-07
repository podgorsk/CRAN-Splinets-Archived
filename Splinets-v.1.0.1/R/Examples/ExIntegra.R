#------------------------------------#
#--- Generate indefinite integral ---#
#------------------------------------#
n=18; k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
# generate a random matrix S
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) #constructing a spline
plot(spl)

dspl = deriva(spl) #derivative
plot(dspl)
is.splinets(dspl)
dintegra(dspl) #the definite integral is 0 (the boundary conditions for 'spl')

ispl = integra(spl) #the integral of a spline
plot(ispl) #the boundary condition on the rhs not satisfied (non-zero value)
ispl@smorder
is.splinets(ispl) #the object does not satisfy the boundary condition for the spline

spll = integra(dspl)
plot(spll)
is.splinets(spll) #the boundary conditions of the integral of the derivative satisfied.


#----------------------------------------------#
#--- Examples with different support ranges ---#
#----------------------------------------------#
n=25; k=2; 
set.seed(5)
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
#Defining support ranges for three splines
supp=matrix(c(2,12,4,20,6,25),byrow=TRUE,ncol=2)
#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[1,2]-supp[1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[2,2]-supp[2,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[3,2]-supp[3,1]+1)*(k+1)),ncol=(k+1)) 
spl=construct(xi,k,SS1,supp[1,]) #constructing the first proper spline
nspl=construct(xi,k,SS2,supp[2,],'CRFC')
spl=gather(spl,nspl) #the second and the first together
nspl=construct(xi,k,SS3,supp[3,],'CRLC')
spl=gather(spl,nspl) #the third is added

plot(spl)
spl@supp


dspl = deriva(spl) #derivative of the splines
plot(dspl)

dintegra(dspl) #the definite integral over the entire range of knots is zero

idspl = integra(dspl) #integral of the derivative returns the original splines
plot(idspl)
is.splinets(idspl)  #and confirms that the object is a spline with boundary conditions 
                    #satified


idspl@supp #Since integral is taken over a function that integrates to zero over 
spl@supp  #each of the support interval, the support of all three objects are the same. 
dspl@supp

ispl=integra(spl)
plot(ispl) #the zero boundary condition at the RHS-end for the splines are not satisfied.
is.splinets(ispl) #thus the object is reported as a non-spline

plot(deriva(ispl))
displ=deriva(ispl) 
displ@supp #Comparison of the supports
spl@supp   
           #Here the integrals have extended support as it is taken from a function
ispl@supp  #that does not integrate to zero.



#---------------------------------------#
#---Example with complicated supports---#
#---------------------------------------#
n=40; xi=seq(0,1,by=1/(n+1)); k=3; 
support=list(matrix(c(2,12,15,27,30,40),ncol=2,byrow = TRUE))
sp=new("Splinets",knots=xi,smorder=k,supp=support) 
m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random
sp1 = is.splinets(sp)[[2]] #Comparison of the corrected and the original 'der' matrices

support=list(matrix(c(2,13,17,30),ncol=2,byrow = TRUE))
sp=new("Splinets",knots=xi,smorder=k,supp=support) 
m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random
sp2 = is.splinets(sp)[[2]] 
sp = gather(sp1,sp2) #a group of two splines
plot(sp)

dsp = deriva(sp) #derivative
plot(dsp)

spl = integra(dsp)
plot(spl) #the spline retrieved
spl@supp  #the supports are retrieved as well
sp@supp
is.splinets(spl) #the proper splinet object that satisfies the boundaries

ispl = integra(sp)
plot(ispl) 
ispl@supp #full support shown by empty list in SLOT 'supp'
is.splinets(ispl) #diagnostic confirms no zeros at the boundaries

spll = deriva(ispl)
plot(spll)
spll@supp

