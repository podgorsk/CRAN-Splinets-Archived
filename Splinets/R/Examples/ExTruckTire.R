#-----------------------------------------------------#
#----------- Plotting the trucktire data -------------#
#-----------------------------------------------------#

#Activating data:
 data(tire)
 data(truck)
 
 matplot(tire[,1],tire[,2:11],type='l',lty=1) #ploting the first 10 tire responses
 
 matplot(truck[,1],truck[,2:11],type='l',lty=1) #ploting the first 10 truck responses
 
 #Projecting truck data into splinet bases
 knots1=seq(0,50, by=2)
 Subtruck= truck[2048:3080,] # selecting the truck data that in the interval[0,50]
 TruckProj=project(as.matrix(Subtruck),knots1)
 
 MeanTruck=matrix(colMeans(TruckProj$coeff),ncol=dim(TruckProj$coeff)[2])
 MeanTruckSp=lincomb(TruckProj$basis,MeanTruck)
 
 plot(MeanTruckSp) #the mean spline of the projections
 
 plot(TruckProj$sp,sID=1:10) #the first ten projections of the functional data
 
 Sigma=cov(TruckProj$coeff)
 Spect=eigen(Sigma,symmetric = TRUE)
 
 plot(Spect$values, type ='l',col='blue', lwd=4 ) #the eigenvalues
 
 EigenTruckSp=lincomb(TruckProj$basis,t(Spect$vec))
 plot(EigenTruckSp,sID=1:5) #the first five largest eigenfunctions
 
 