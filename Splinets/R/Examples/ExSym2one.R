#-----------------------------------------------------#
#-------Representations of derivatives at knots-------#
#-----------------------------------------------------#
n=10; k=3; xi=seq(0,1,by=1/(n+1)) #the even number of equally spaced knots 
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) #construction of a spline
a=spl@der[[1]]
b=sym2one(a)
aa=sym2one(b,inv=TRUE) # matrix 'aa' is the same as 'a'

n=11; xi2=seq(0,1,by=1/(n+1)) #the odd number of knots case
S2=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl2=construct(xi2,k,S2) #construction of a spline
a2=spl2@der[[1]]
b2=sym2one(a2)
aa2=sym2one(b2, inv=TRUE) # matrix 'aa2' is the same as 'a2'

#-----------------------------------------------------#
#--------------More complex support sets--------------#
#-----------------------------------------------------#
#Zero order splines, non-equidistant case, support with three components
n=43; xi=seq(0,1,by=1/(n+1)); k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1;
support=list(matrix(c(2,14,17,30,32,43),ncol=2,byrow = TRUE))
#Third order splines
ssp=new("Splinets",knots=xi,supp=support,smorder=k) #with partial support

m=sum(ssp@supp[[1]][,2]-ssp@supp[[1]][,1]+1) #the total number of knots in the support
ssp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1)))  #the derivative matrix at random
IS=is.splinets(ssp) 
IS$robject@der
IS$robject@supp
b=sym2one(IS$robject@der[[1]],IS$robject@supp[[1]]) #the RHS limits at the knots
a=sym2one(b,IS$robject@supp[[1]],inv=TRUE) #is the same as the SLOT supp in IS@robject