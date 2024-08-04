#--------------------------------------#
#----Splinet, equally spaced knots-----#
#--------------------------------------#
k=2 # order
n_knots = 12 # number of knots
xi = seq(0, 1, length.out = n_knots)

so = splinet(xi, k)

plot(so$bs) #Plotting B-splines
plot(so$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(so$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#An example of the dyadic structure with equally spaced knots
k=3
N=3
n_knots=2^N*k-1 #the number of internal knots for the dyadic case

xi = seq(0, 1, length.out = n_knots+2)

so = splinet(xi) 

plot(so$bs,type="simple",vknots=FALSE,lwd=3) #Plotting B-splines in a single simple plot
plot(so$os,type="simple",vknots=FALSE,lwd=3) 

plot(so$os,lwd=3,mrgn=2) #Plotting the splinet on the dyadic net of support intervals

so=splinet(xi, Bsplines=so$bs, type='gsob') #Obtaining the Gram-Schmidt orthogonalization
plot(so$os,type="simple",vknots=FALSE)      #Without computing B-splines again

so=splinet(xi, Bsplines=so$bs, type='twob') #Obtaining the symmetrize orthogonalization
plot(so$os,type="simple",vknots=FALSE)     

#-------------------------------------#
#---Splinet, unequally spaced knots---#
#-------------------------------------#
n_knots=30

xi = c(0, sort(runif(n_knots)), 1)

sone = splinet(xi, k)

plot(sone$bs, type='dyadic') #Plotting B-splines
plot(sone$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sone$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#------------------------------------------#
#---Dyadic splinet, equally spaced knots---#
#------------------------------------------#
k = 2 # order
N = 3 # support level
n_so = k*(2^N-1) # number of splines in a dyadic structure with N and k
n_knots = n_so + k + 1 # number of knots
xi = seq(0, 1, length.out = n_knots)

sodyeq = splinet(xi, k)

plot(sodyeq$bs) #Plotting B-splines
plot(sodyeq$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sodyeq$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#--------------------------------------------#
#---Dyadic splinet, unequally spaced knots---#
#--------------------------------------------#
xi = c(0, sort(runif(n_knots)), 1)

sody = splinet(xi, k)


plot(sody$bs) #Plotting B-splines
plot(sody$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sody$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#-----------------------------------------#
#---Bspline basis, equally spaced knots---#
#-----------------------------------------#
n = 15
xi = seq(0,1,length.out = n+2)
order = 2

bs = splinet(xi, order, type = 'bs')

plot(bs$bs)

#---------------------------------------------#
#---Bspline basis, non-equally spaced knots---#
#---------------------------------------------#
n = 6
xi = c(0,sort(runif(n)),1)
order = 3

so = splinet(xi, order, type = 'bs') #unnormalized version
plot(so$bs) 

so1 = splinet(type='bs',Bsplines=so$bs,norm=TRUE) #normalized version
plot(so1$bs)

#-------------------------------------------------#
#---Gram-Schmidt osplines, equally spaced knots---#
#-------------------------------------------------#
so = splinet(xi, order,  type = 'gsob')

plot(so$bs)
plot(so$os)

#Using the previously generated B-splines and normalizing them
so1 = splinet(Bsplines=so$bs, type = "gsob",norm=TRUE) 

plot(so1$bs) #normalized B-splines
plot(so1$os) #the one sided osplines

gm = gramian(so1$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm))) #verification of the orghonoalization of the matrix

#-----------------------------------------------------#
#---Gram-Schmidt osplines, non-equally spaced knots---#
#-----------------------------------------------------#

so = splinet(Bsplines=sody$bs, type = 'gsob') #previously genereted Bsplines

plot(so$bs)
plot(so$os)

gm = gramian(so$os)
diag(gm)
sum(gm - diag(diag(gm)))

#---------------------------------------------#
#---Twosided osplines, equally spaced knots---#
#---------------------------------------------#

so = splinet(Bsplines=bs$bs, type = 'twob')
plot(so$os)

gm = gramian(so$os) #verification of the orthogonality
diag(gm)
sum(gm - diag(diag(gm)))

#-------------------------------------------------#
#---Twosided osplines, non equally spaced knots---#
#-------------------------------------------------#

so = splinet(Bsplines=sody$bs, type = 'twob')
plot(so$os)

gm = gramian(so$os) #verification of the orthogonality
diag(gm)
sum(gm - diag(diag(gm)))


