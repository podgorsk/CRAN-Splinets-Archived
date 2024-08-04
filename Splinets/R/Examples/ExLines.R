#-----------------------------------------------------#
#------Adding spline lines to an existing graph-------#
#-----------------------------------------------------#
n=17; k=4; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) 
plot(spl,main="Mean Spline",lty=2,lwd=2)

RS=rspline(spl,5)
plot(RS,main="Random splines around the mean spline" )
lines(spl,col='red',lwd=4,lty=2)

