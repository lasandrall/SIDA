sidanettunerange=function(Xdata=Xdata,Y=Y,ngrid=8,standardize=TRUE,weight=0.5,eta=0.5,myedges=myedges,myedgeweight=myedgeweight,withCov=FALSE){


#check size of each data

dsizes=lapply(Xdata, function(x) dim(x))
n=dsizes[[1]][1]
p=lapply(Xdata, function(x) dim(x)[2])
D=length(dsizes)

if(is.null(withCov)){
   withCov=FALSE
}

if(is.null(ngrid)){
   ngrid=8
}

if(is.null(standardize)){
   standardize=TRUE
}

if(is.null(weight)){
   weight=0.5
}

if(is.null(eta)){
   eta=1e-03
}

if(withCov==TRUE){
  D=D-1
}

nK=length(unique(as.vector(Y))) -1

#standardize if true
if(standardize==TRUE){
   Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
}

#create normalized Laplacian of Graph
mynormLaplacianG=myNLaplacianG(Xdata,myedges,myedgeweight)
#obtain nonsparse solutions
mynsparse=myfastIDAnonsparse(Xdata,Y,weight)
myfinner=myfastinner(Xdata,Y,mynsparse$sqrtminvmat,mynsparse$myalphaoldmat,mynsparse$tildealphamat, weight)


#obtain upper and lower bounds
ubx=lapply(myfinner$SepAndAssocd, function(x) norm(x,'i')/1.2)
lbx=lapply(1:D, function(x) 1.2*sqrt(log(p[[x]])/n)*ubx[[x]])
ubx=lapply(1:D, function(x) ubx[[x]])

#tuning range for each data
Taugrid=list()
cc=lapply(1, function(x1,x2)  cbind(lbx,ubx))

cc=as.matrix(do.call(rbind,cc))
for(d in 1:D){
  Taugrid[[d]]=seq(as.numeric(cc[d,1]),as.numeric(cc[d,2]),length.out=(ngrid+1))
}

#about 25% sparsity

myperx=lapply(Taugrid, function(x) quantile(x[1:ngrid], c(.1, .15, .2, .25, .35, .45), type=5))#similar to matlab
myperx2=do.call(rbind,myperx)
 for(loc in 1:6){
   mTaux=sapply(1:D, function(i) list(t(myperx2[i,loc])))
   myres=sidanetinner(Xdata,Y,mynsparse$sqrtminvmat,mynsparse$myalphaoldmat,mynsparse$tildealphamat, mynsparse$tildelambda,mTaux,weight,eta,myedges,myedgeweight,mynormLaplacianG,withCov)
   nnz=sapply(1:D, function(i) list(colSums(myres$hatalpha[[i]]!=0)/dsizes[[i]][2]))
   nnz=cbind(c(do.call(rbind,nnz))) #turns into a vector
   if(all(nnz<=0.25)){
     break
   }
}
#final grid
Tauvec=sapply(1:D, function(i) list(seq(as.numeric(t(myperx2[i,loc])),as.numeric(ubx[[i]]),len=(ngrid+1))))
Tauvec=sapply(1:D, function(x) list(Tauvec[[x]][1:ngrid]))

result=list(Tauvec=Tauvec)
return(result)
}

