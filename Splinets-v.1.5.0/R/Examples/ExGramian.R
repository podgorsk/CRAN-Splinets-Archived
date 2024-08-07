#---------------------------------------#
#---- Simple three splines example -----# 
#---------------------------------------#
n=25; k=3
xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
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

plot(spl)
gramian(spl)
gramian(spl, norm_only = TRUE)
gramian(spl, sID = c(1,3))
gramian(spl,sID=c(2,3),Sp2=spl,s2ID=c(1)) #the cross-Gramian matrix  

#-----------------------------------------#
#--- Example with varying support sets ---#
#-----------------------------------------#
n=40; xi=seq(0,1,by=1/(n+1)); k=2; 
support=list(matrix(c(2,9,15,24,30,37),ncol=2,byrow = TRUE))
sp=new("Splinets",knots=xi,smorder=k,supp=support) 
m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random
sp1 = is.splinets(sp)[[2]] #the correction of 'der' matrices

support=list(matrix(c(5,12,17,29),ncol=2,byrow = TRUE))
sp=new("Splinets",knots=xi,smorder=k,supp=support) 
m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random
sp2 = is.splinets(sp)[[2]] 

spp = gather(sp1,sp2)

support=list(matrix(c(3,10,14,21,27,34),ncol=2,byrow = TRUE))
sp=new("Splinets",knots=xi,smorder=k,supp=support) 
m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random
sp3 = is.splinets(sp)[[2]] 

spp = gather(spp, sp3)

plot(spp)
gramian(spp) #the regular gramian matrix
spp2=subsample(spp,sample(1:3,size=3,rep=TRUE))
gramian(Sp=spp,Sp2=spp2) #cross-Gramian matrix


#-----------------------------------------#
#--------- Grammian for B-splines --------#
#-----------------------------------------#

n=25; xi=seq(0,1,by=1/(n+1)); k=2; 
Sp=splinet(xi) #B-splines and corresponding splinet

gramian(Sp$bs) #band grammian matrix for B-splines
gramian(Sp$os) #diagonal gramian matrix for the splinet
A=gramian(Sp=Sp$bs,Sp2=Sp$os) #cross-Gramian matrix, the coefficients of
                              #the decomposition of the B-splines

plot(Sp$bs)
plot(lincomb(Sp$os,A))