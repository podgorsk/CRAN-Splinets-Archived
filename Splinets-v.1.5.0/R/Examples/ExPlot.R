#-----------------------------------------------------#
#-------------------Ploting splinets------------------#
#-----------------------------------------------------#
#Constructed splines

n=25; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1; k=3
supp=list(t(c(2,12)),t(c(4,20)),t(c(6,25))) #defining support ranges for three splines

#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[[1]][1,2]-supp[[1]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[[2]][1,2]-supp[[2]][1,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[[3]][1,2]-supp[[3]][1,1]+1)*(k+1)),ncol=(k+1)) 

spl=construct(xi,k,SS1,supp[[1]]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[[2]],'CRFC')
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[[3]],'CRLC')
spl=gather(spl,nspl) #the third is added

plot(spl)
  plot(spl,sID=c(1,3))
plot(spl,sID=2)
t = seq(0,0.5,length.out = 1000)
plot(spl, t, sID = 1)

#Random splines
n=17; k=4; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) 
plot(spl,main="Mean Spline",lty=2,lwd=2,xlab='')

RS=rspline(spl,5)
plot(RS,main="Random splines around the mean spline",ylim=3*range(spl@der[[1]][,1]) )
lines(spl,col='red',lwd=4,lty=2)

#Periodic splines
xi = seq(0, 1, length.out = 25)
so = splinet(xi, periodic = TRUE) 

plot(so$bs)
plot(so$os)
plot(so$bs,type= "dyadic")
plot(so$bs, sID=c(4,6))
plot(so$os, type="simple",sID=c(4,6))


