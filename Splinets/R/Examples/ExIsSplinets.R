#-----------------------------------------------------#
#--------Diagnostics of simple Splinets objects-------#
#-----------------------------------------------------#
#----------Full support equidistant cases-------------#
#-----------------------------------------------------#
#Zero order splines, equidistant case, full support
n=20; xi=seq(0,1,by=1/(n+1))
sp=new("Splinets",knots=xi) 
sp@equid #equidistance flag
#Diagnostic of 'Splinets' object 'sp'
is.splinets(sp)

IS=is.splinets(sp)
IS[[1]] #informs if the object is a spline
IS$is   #equivalent to the above

#Third order splines with a noisy matrix of the derivative
set.seed(5)
k=3; sp@smorder=k; sp@der[[1]]=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))

IS=is.splinets(sp)
IS[[2]]@taylor #corrections
sp@taylor
IS[[2]]@der #corrections
sp@der

is.splinets(IS[[2]]) #The output object is a valid splinet

#-----------------------------------------------------#
#--------Full support non-equidistant cases-----------#
#-----------------------------------------------------#
#Zero order splines, non-equidistant case, full support
set.seed(5)
n=17; xi=sort(runif(n+2))
xi[1]=0 ;xi[n+1]=1 #The last knot is not in the order.
          #It will be reported and corrected in the output.  
sp=new("Splinets",knots=xi) 

xi #original knots
sp@knots #vs. corrected ones
sp@taylor

#Diagnostic of 'Splinets' object 'sp'
is.splinets(sp)

IS=is.splinets(sp)
 
nsp=IS$robject #the output spline -- a corrected version of the input
nsp@der
sp@der

#Third order splines
nsp@smorder=3
IS=is.splinets(nsp)

IS[[2]]@taylor #corrections
nsp@taylor
IS[[2]]@der #corrections
nsp@der
is.splinets(IS[[2]]) #verification that the correction is a valid object

#Randomly assigning the derivative -- a very 'unstable' function.
set.seed(5)
k=nsp@smorder; S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1)); nsp@der[[1]]=S

IS=is.splinets(nsp) #the 2nd element of 'IS' is a spline obtained by correcting 'S'
nsp=is.splinets(IS[[2]])
nsp$is #The 'Splinets' object is correct, alternatively use 'nsp$[[1]]'.
nsp$robject #A correct spline object, alternatively use 'nsp$[[2]]'.

#-----------------------------------------------------#
#-----Splinets objects with varying support sets------#
#-----------------------------------------------------#
#-----------------Eequidistant cases------------------#
#-----------------------------------------------------#
#Zero order splines, equidistant case, support with three components
n=20; xi=seq(0,1,by=1/(n+1))
support=list(matrix(c(2,5,6,8,12,18),ncol=2,byrow = TRUE)) 
sp=new("Splinets",knots=xi,supp=support) 
is.splinets(sp)

IS=is.splinets(sp)

sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
dim(IS[[2]]@der[[1]])[1]  #the number of rows in the derivative matrix
IS[[2]]@der[[1]] #the corrected object
sp@der #the input derivative matrix

#Third order splines
n=40; xi=seq(0,1,by=1/(n+1)); k=3; 
support=list(matrix(c(2,12,15,27,30,40),ncol=2,byrow = TRUE))
m=sum(support[[1]][,2]-support[[1]][,1]+1) #the number of knots in the support
SS=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random

sp=new("Splinets",knots=xi,smorder=k,supp=support,der=SS) 

IS=is.splinets(sp)

m=sum(sp@supp[[1]][,2]-sp@supp[[1]][,1]+1) #the number of knots in the support
sp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1))) #the derivative matrix at random

IS=is.splinets(sp) #Comparison of the corrected and the original 'der' matrices
sp@der
IS[[2]]@der
is.splinets(IS[[2]]) #verification

#-----------------------------------------------------#
#----------------Non-equidistant cases----------------#
#-----------------------------------------------------#
#Zero order splines, non-equidistant case, support with three components
set.seed(5)
n=43; xi=seq(0,1,by=1/(n+1)); k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1;
support=list(matrix(c(2,14,17,30,32,43),ncol=2,byrow = TRUE))

ssp=new("Splinets",knots=xi,supp=support) #with partial support
nssp=is.splinets(ssp)$robject
nssp@supp
nssp@der

#Third order splines
nssp@smorder=3 #changing the order of the 'Splinets' object
set.seed(5)
m=sum(nssp@supp[[1]][,2]-nssp@supp[[1]][,1]+1) #the number of knots in the support
nssp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1)))  #the derivative matrix at random
IS=is.splinets(nssp) 
IS$robject@der
is.splinets(IS$robject)$is #verification of the corrected output object

