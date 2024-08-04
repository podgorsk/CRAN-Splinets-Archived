#-------------------------------------------------------#
#--The support layers of the dyadic structure of bases--#
#-------------------------------------------------------#
k=4 # order
n = 36 # number of the internal knots (excluding the endpoints)
xi = seq(0, 1, length.out = n+2)
spnt=splinet(xi,k)

plot(spnt$os)                #standard plotting 
plot(spnt$bs,type='dyadic')  #dyadic format of plots

net=seq2dyad(n-k+1,k)  #retrieving dyadic structure
ind1=c(net[[4]][[1]],net[[4]][[2]])

plot(subsample(spnt$os,ind1))

ind2=c(net[[4]][[3]],net[[4]][[4]]) #the lowest support in the dyadic net

lines(subsample(spnt$bs,ind2))