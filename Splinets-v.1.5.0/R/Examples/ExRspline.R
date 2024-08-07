#-----------------------------------------------------#
#-------Simulation of a standard random splinet-------#
#-----------------------------------------------------#
n=17; k=4; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) #Construction of the mean spline

RS=rspline(spl)
graphsp=evspline(RS) #Evaluating the random spline
meansp=evspline(spl)
RS=rspline(spl,5) #Five more samples
graphsp5=evspline(RS)

m=min(graphsp[,2],meansp[,2],graphsp5[,2:6])
M=max(graphsp[,2],meansp[,2],graphsp5[,2:6])

plot(graphsp,type='l',ylim=c(m,M))
lines(meansp,col='red',lwd=3,lty=2) #the mean spline
for(i in 1:5){lines(graphsp5[,1],graphsp5[,i+1],col=i)}

#-----------------------------------------------------#
#------------Different construction method------------#
#-----------------------------------------------------#
RS=rspline(spl,8,mthd='CRLC'); graphsp8=evspline(RS)

m=min(graphsp[,2],meansp[,2],graphsp8[,2:6])
M=max(graphsp[,2],meansp[,2],graphsp8[,2:6])

plot(meansp,col='red',type='l',lwd=3,lty=2,ylim=c(m,M)) #the mean spline
for(i in 1:8){lines(graphsp8[,1],graphsp8[,i+1],col=i)}

#-----------------------------------------------------#
#-------Simulation of with different variances--------#
#-----------------------------------------------------#
Sigma=seq(0.1,1,n+2);Theta=seq(0.1,1,k+1)
RS=rspline(spl,N=10,Sigma=Sigma) #Ten samples
RS2=rspline(spl,N=10,Sigma=Sigma,Theta=Theta) #Ten samples
graphsp10=evspline(RS); graphsp102=evspline(RS2)

m=min(graphsp[,2],meansp[,2],graphsp10[,2:10])
M=max(graphsp[,2],meansp[,2],graphsp10[,2:10])

plot(meansp,type='l',ylim=c(m,M),col='red',lwd=3,lty=2) 
for(i in 1:10){lines(graphsp10[,1],graphsp10[,i+1],col=i)}

m=min(graphsp[,2],meansp[,2],graphsp102[,2:10])
M=max(graphsp[,2],meansp[,2],graphsp102[,2:10])

plot(meansp,type='l',ylim=c(m,M),col='red',lwd=3,lty=2) 
for(i in 1:10){lines(graphsp102[,1],graphsp102[,i+1],col=i)}

#-----------------------------------------------------#
#-------Simulation for the mean spline to be----------#
#------=----defined on incomplete supports------------#
#-----------------------------------------------------#
n=43; xi=seq(0,1,by=1/(n+1)); k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1;
support=list(matrix(c(2,14,25,43),ncol=2,byrow = TRUE))
ssp=new("Splinets",knots=xi,supp=support) #with partial support
nssp=is.splinets(ssp)$robject
nssp@smorder=3 #changing the order of the 'Splinets' object
m=sum(nssp@supp[[1]][,2]-nssp@supp[[1]][,1]+1) #the number of knots in the support
nssp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1)))  #the derivative matrix at random
spl=is.splinets(nssp)$robject
RS=rspline(spl,Sigma=0.05,Theta=c(1,0.5,0.3,0.05))
graphsp=evspline(RS);
meansp=evspline(spl)

m=min(graphsp[,2],meansp[,2],graphsp5[,2:6])
M=max(graphsp[,2],meansp[,2],graphsp5[,2:6])

plot(graphsp,type='l',ylim=c(m,M))
lines(meansp,col='red',lwd=3,lty=2) #the mean spline

