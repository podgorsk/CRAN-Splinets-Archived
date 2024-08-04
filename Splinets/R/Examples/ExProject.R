#-------------------------------------------------#
#----Representing splines in the spline bases-----#
#-------------------------------------------------#
k=3 # order
n = 10 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 

plot(spl) # plotting a spline
spls=rspline(spl,5) # a random sample of splines

Repr=project(spl) #decomposition of splines into the splinet coefficients

Repr=project(spl, graph = TRUE) #decomposition of splines with the following graphs
                             #that illustrate the decomposition: 
                             # 1) The orthogonal spine basis on the dyadic grid; 
                             # 2) The coefficients of the projections on the dyadic grid;
                             # 3) The input splines;
                             # 4) The projections of the input.
                
Repr$coeff       #the coefficients of the decomposition
plot(Repr$sp) #plot of the reconstruction of the spline

plot(spls)
Reprs=project(spls,basis = Repr$basis) #decomposing splines using the available basis
plot(Reprs$sp) 

Reprs2=project(spls,type = 'gsob') #using the Gram-Schmidt basis


#The case of the regular non-normalized B-splines:
Reprs3=project(spls,type = 'bs') 
plot(Reprs3$basis)
gramian(Reprs3$basis,norm_only = TRUE) #the B-splines follow the classical definition and 
                                       #thus are not normalized
plot(spls)

plot(Reprs3$basis) #Bsplines
plot(Reprs3$sp) #reconstruction using the B-splines and the decomposition

#a non-equidistant example
n=10; k=3
set.seed(5)
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1 
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) 
plot(spl)
spls=rspline(spl,5) # a random sample of splines
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

spls=rspline(spl,5) # a random sample of splines
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

k=3 # order
n = 10 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

spl=construct(xi,k,S) 
spls=rspline(spl,5) # a random sample of splines

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
k=3; n = 10; xi = seq(0, 1, length.out = n+2)
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S); spls=rspline(spl,10) # a random sample of splines

x=runif(50)
FData=evspline(spls,x=x) #discrete functional data

matplot(FData[,1],FData[,-1],pch='.',cex=3)
                   #adding small noise to the data
noise=matrix(rnorm(length(x)*10,0,sqrt(var(FData[,2]/10))),ncol=10)

FData[,-1]=FData[,-1]+noise

matplot(FData[,1],FData[,-1],pch='.',cex=3)

knots=runif(12) 

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

#Using two-sided orthonormal basis

DatProj5=project(FData,knots,type='twob')
matplot(FData[,1],FData[,-1],pch='.',cex=3)
lines(DatProj5$sp)
lines(spls)

#--------------------------------------------------#
#-----Projecting into a periodic spline space------#
#--------------------------------------------------#
#generating periodic splines
n=1# number of samples
k=3
N=3

n_knots=2^N*k-1 #the number of internal knots for the dyadic case
xi = seq(0, 1, length.out = n_knots+2)

so = splinet(xi,smorder = k, periodic = TRUE) #The splinet basis
stwo = splinet(xi,smorder = k,type='twob', periodic = TRUE) #The two-sided orthogonal basis

plot(so$bs,type='dyadic',main='B-Splines on dyadic structure') #B-splines on the dyadic graph 

plot(stwo$os,main='Symmetric OB-Splines') #The two-sided orthogonal basis 

plot(stwo$os,type='dyadic',main='Symmetric OB-Splines on dyadic structure')

# generating a periodic spline as a linear combination of the periodic splines 
A1= as.matrix(c(1,0,0,0.7,0,0,0,0.8,0,0,0,0.4,0,0,0, 1, 0,0,0,0,0,1,0, .5),nrow= 1)
circular_spline=lincomb(so$os,t(A1))

plot(circular_spline)

#Graphical visualizations of the projections 

Pro_spline=project(circular_spline,basis = so$os,graph = TRUE)
plot(Pro_spline$sp)

#---------------------------------------------------------------#
#---Projecting discretized data into a periodic spline space----#
#---------------------------------------------------------------#
nx=100 # number of discritization 
n=1# number of samples
k=3
N=3

n_knots=2^N*k-1 #the number of internal knots for the dyadic case
xi = seq(0, 1, length.out = n_knots+2)

so = splinet(xi,smorder = k, periodic = TRUE)

hf=1/nx 
grid=seq (hf , 1, by=hf) #grid 
l=length(grid)
BB = evspline(so$os, x =grid)
fbases=matrix(c(BB[,2],BB[,5],BB[,9],BB[,13],BB[,17], BB[,23], BB[,25]), nrow = nx)

#constructing periodic data 
f_circular=matrix(0,ncol=n+1,nrow=nx)
lambda=c(1,0.7,0.8,0.4, 1,1,.5)
f_circular[,1]= BB[,1]
f_circular[,2]= fbases%*%lambda

plot(f_circular[,1], f_circular[,2], type='l')

Pro=project(f_circular,basis = so$os) 
plot(Pro$sp)



