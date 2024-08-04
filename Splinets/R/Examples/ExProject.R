#-------------------------------------------------#
#----Representing splines in the spline bases-----#
#-------------------------------------------------#
\donttest{
k=3 # order
n = 16 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 

plot(spl) # plotting a spline
spls=rspline(spl,10) # a random sample of splines

Repr=project(spl) #decomposition of splines into the splinet coefficients

plot(Repr$basis) #the plot of the splinet used for decomposition
Repr$coeff       #the coefficients of the decomposition
plot(Repr$sp) #plot of the reconstruction of the spline

plot(spls)
Reprs=project(spls,basis = Repr$basis) #decomposing splines using the available basis
plot(Reprs$sp) 

Reprs2=project(spls,type = 'gsob') #using the Gram-Schmidt basis

plot(Reprs2$basis)
plot(Reprs2$sp) 

#The case of the regular non-normalized B-splines:
Reprs3=project(spls,type = 'bs') 
plot(Reprs3$basis)
gramian(Reprs3$basis,norm_only = TRUE) #the B-splines follow the classical definition and 
                                       #thus are not normalized
plot(spls)

plot(Reprs3$basis) #Bsplines
plot(Reprs3$sp) #reconstruction using the B-splines and the decomposition

#a non-equidistant example
n=17; k=4
set.seed(5)
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) 
plot(spl)
spls=rspline(spl,20) # a random sample of splines
plot(spls)
Reprs=project(spls,type = 'twob') #decomposing using the two-sided orthogonalization
plot(Reprs$basis)
plot(Reprs$sp) 

#The case of the regular non-normalized B-splines:
Reprs2=project(spls,basis=Reprs$basis) 
plot(Reprs2$sp) #reconstruction using the B-splines and the decomposition

#-------------------------------------------------#
#-----Projecting splines into a spline space------#
#-------------------------------------------------#

k=3 # order
n = 10 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 

plot(spl) #the spline

knots=runif(8)
Prspl=project(spl,knots)

plot(Prspl$sp) #the projection spline
Rspl=refine(spl,newknots = knots) #embedding the spline to the common space
plot(Rspl)
RPspl=refine(Prspl$sp,newknots = xi)  #embedding the projection spline to the common space
plot(RPspl)
All=gather(RPspl,Rspl) #creating the Splinets-object with the spline and its projection
Rbasis=refine(Prspl$basis,newknots = xi) #embedding the basis to the common space
plot(Rbasis)

Res=lincomb(All,matrix(c(1,-1),ncol=2))
plot(Res)
gramian(Res,Sp2 = Rbasis) #the zero valued innerproducts -- the orthogonality of the residual spline

spls=rspline(spl,10) # a random sample of splines
Prspls=project(spls,knots,type='bs') #projection in the B-spline representation
plot(spls)
lines(Prspls$sp) #presenting projections on the common plot with the original splines

Prspls$sp@knots
Prspls$sp@supp

plot(Prspls$basis)   #Bspline basis


#An example with partial support

Bases=splinet(xi,k)

BS_Two=subsample(Bases$bs,c(2,length(Bases$bs@der)))
plot(BS_Two)
A=matrix(c(1,-2),ncol=2)
spl=lincomb(BS_Two,A)
plot(spl)

knots=runif(13)
Prspl=project(spl,knots)
plot(Prspl$sp)
Prspl$sp@knots
Prspl$sp@supp

#Using explicit bases 

k=5 # order
n = 40 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 
spls=rspline(spl,20) # a random sample of splines

plot(spls)

knots=runif(20)
base=splinet(knots,smorder=k)
plot(base$os)

Prsps=project(spls,basis=base$os)
plot(Prsps$sp) #projection splines vs. the original splines
lines(spls)

#------------------------------------------------------#
#---Projecting discretized data into a spline space----#
#------------------------------------------------------#

k=4; n = 30; xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S); spls=rspline(spl,10) # a random sample of splines

x=runif(80)
FData=evspline(spls,x=x) #discrete functional data

matplot(FData[,1],FData[,-1],pch='.',cex=3)
                   #adding small noise to the data
noise=matrix(rnorm(length(x)*10,0,sqrt(var(FData[,2]/10))),ncol=10)

FData[,-1]=FData[,-1]+noise

matplot(FData[,1],FData[,-1],pch='.',cex=3)

knots=runif(14) 

DatProj=project(FData,knots)

lines(DatProj$sp) #the projections at the top of the original noised data

plot(DatProj$basis) #the splinet in the problem

#Adding knots to the projection space so that all data points are included
#in the range of the knots for the splinet basis

knots=c(-0.1,0.0,0.1,0.85, 0.9, 1.1,knots)

bases=splinet(knots)

DatProj1=project(FData,basis = bases$os)


matplot(FData[,1],FData[,-1],pch='.',cex=3)
lines(DatProj1$sp) 

#Using the B-spline basis

knots=xi

bases=splinet(knots,type='bs')

DatProj3=project(FData,basis = bases$bs)

matplot(FData[,1],FData[,-1],pch='.',cex=3)
lines(DatProj3$sp) 

DatProj4=project(FData,knots,k,type='bs') #this includes building the base of order 4

matplot(FData[,1],FData[,-1],pch='.',cex=3)
lines(DatProj4$sp) 
lines(spls) #overlying the functions that the original data were built from

plot(DatProj3$basis) #the 3rd order
plot(DatProj4$basis) #the 4th order

#Using two-sided orthonormal basis

DatProj5=project(FData,knots,type='twob')
matplot(FData[,1],FData[,-1],pch='.',cex=3)
lines(DatProj5$sp)
lines(spls)

plot(DatProj5$basis) #two-sided orthonormal spline basis


#-------------------------------------------------#
#---------Simple functional data analysis---------#
#-------------------------------------------------#

matplot(truck[,1],truck[,2:dim(truck)[2]],type='l',lty=1)

knots=seq(-100,100, by=1)
TruckProj=project(as.matrix(truck),knots)

MeanTruck=matrix(colMeans(TruckProj$coeff),ncol=dim(TruckProj$coeff)[2])
MeanTruckSp=lincomb(TruckProj$basis,MeanTruck)

plot(MeanTruckSp) #the mean spline of the projections

plot(TruckProj$sp,sID=1:10) #the first ten projections of the functional data

Sigma=cov(TruckProj$coeff)
Spect=eigen(Sigma,symmetric = TRUE)

plot(Spect$values, type ='l',col='blue', lwd=4 ) #the eigenvalues

EigenTruckSp=lincomb(TruckProj$basis,t(Spect$vec))
plot(EigenTruckSp,sID=1:5) #the first five largest eigenfunctions

#-------------------------------------------------#
#------------Graphical features ------------------#
#-------------------------------------------------#

project(FData,knots, graph = TRUE) #the splinet projection

project(FData,basis = bases$bs,graph = TRUE) #B-spline projection

project(FData,knots,type='twob',graph = TRUE) #two-sided projection
}
