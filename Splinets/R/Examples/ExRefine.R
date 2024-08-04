#-------------------------------------------------#
#----Refining splines - the full support case-----#
#-------------------------------------------------#
k=3 # order
n = 16 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 

plot(spl) # plotting a spline


rspl=refine(spl) # refining the equidistant by doubling its knots
plot(rspl)
rspl@equid  # the outcome is equidistant

#a non-equidistant case
n=17; k=4
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) 
plot(spl)
mult=3 #adding two knots between each subsequent pair of the original knots
rspl=refine(spl,mult) 
is.splinets(rspl)
plot(rspl)

#adding specific knots
rspl=refine(spl,newknots=c(0.5,0.75))
rspl@knots
is.splinets(rspl)
plot(rspl)

#----------------------------------------------------#
#----Refining splines - the partial support case-----#
#----------------------------------------------------#

Bases=splinet(xi,k)
plot(Bases$bs)
Base=Bases$bs

BS_Two=subsample(Bases$bs,c(1,length(Base@der)))
plot(BS_Two)
A=matrix(c(1,-1),ncol=2)
spl=lincomb(BS_Two,A)

rspl=refine(spl) #doubling the number of knots
plot(rspl)
is.splinets(rspl)
rspl@supp #the support is evaluated 
spl@supp

#The case of adding knots explicitely
BS_Middle=subsample(Bases$bs,c(floor(length(Base@der)/2)))
spls=gather(spl,BS_Middle)
plot(spls)

rspls=refine(spls, newknots=c(0.2,0.5,0.85)) #two splines with partial support sets 
                                             #by adding three knots to B-splines
plot(rspls)

#----------------------------------------------------#
#------Refining splines over the larger range--------#
#----------------------------------------------------#

k=4 # order
n = 25 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 

plot(spl) # plotting a spline
newknots=c(-0.1,0.4,0.6,1.2)  #the added knots create larger range 
rspl=refine(spl,newknots=newknots)
spl@supp #the original spline has the full support
rspl@supp #the embedded spline has partial support
spl@equid
rspl@equid
plot(rspl)
