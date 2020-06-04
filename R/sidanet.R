sidanet=function(Xdata,Y,myedges,myedgeweight,Tau,withCov=FALSE,Xtestdata=NULL,Ytest=NULL,AssignClassMethod='Joint',plotIt=FALSE,standardize=TRUE,maxiteration=20,weight=0.5,thresh= 1e-03,eta=0.5,mynormLaplacianG=NULL){


  #check inputs640
  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  nsizes=lapply(Xdata, function(x) dim(x)[1])

  if(all(nsizes!=nsizes[[1]])){
    stop('The datasets  have different number of observations')
  }

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }

  if(is.null(Xtestdata)){
    Xtestdata=Xdata
  }
  
  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])
  if(all(ntestsizes!=ntestsizes[[1]])){
    stop('The testing datasets  have different number of observations')
  }
  
  if(is.null(Ytest)){
    Ytest=Y
  }
  
  if(is.null(withCov)){
    withCov=FALSE
  }
  
  if(is.null(plotIt)){
    plotIt=FALSE
  }
  
  if(is.null(standardize)){
    standardize=TRUE
  }
  
  if(is.null(maxiteration)){
    maxiteration=20
  }
  
  if(is.null(weight)){
    weight=0.5
  }
  
  if(is.null(thresh)){
    thresh=1e-03
  }
  
  if(is.null(mynormLaplacianG)){
    mynormLaplacianG=myNLaplacianG(Xdata,myedges,myedgeweight)
  }
  
 #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }


  nK=length(unique(as.vector(Y))) -1


  #norm function for convergence
  normdiff=function(xnew,xold){
    ndiff=norm(xnew-xold,'f')^2 / norm (xold,'f')^2
  }
  #initialize
  iter=0
  diffalpha=1
  reldiff=1
  relObj=1

  mynsparse=myfastIDAnonsparse(Xdata,Y,weight)
  myalpha=mynsparse$myalphaoldmat
  #while convergence is not met
  while(iter < maxiteration && min(min(reldiff,relObj),max(diffalpha))> thresh){
    iter=iter+1
    print(iter)
    myalphaold=myalpha
    mysidainner=sidanetinner(Xdata,Y,mynsparse$sqrtminvmat,myalphaold,mynsparse$tildealphamat, mynsparse$tildelambda,Tau,weight,eta,myedges,myedgeweight,mynormLaplacianG,withCov)
    
    myalpha=mysidainner$hatalpha
    nz=sapply(1:D, function(i) list(colSums(myalpha[[i]]!=0)))
    nz=cbind(c(do.call(rbind,nz)))
    if(any(nz==0)){
      myalpha=myalphaold
      break
    }
    diffalpha=mapply(normdiff, myalpha, myalphaold)
    sumnormdiff=sum(sapply(1:D, function(d) norm(myalpha[[d]]-myalphaold[[d]],'f')^2 ))
    sumnormold=sum(sapply(1:D, function(d) norm(myalphaold[[d]],'f')^2    ))
    reldiff=sumnormdiff/sumnormold
    ObjNew=sum(sapply(1:D, function(d) (1-eta)*norm(myalpha[[d]],'f')^2 + eta*norm(as.matrix(mysidainner$myL[[d]])%*%myalpha[[d]],'f')^2))
    ObjOld=sum(sapply(1:D, function(d) (1-eta)*norm(myalphaold[[d]],'f')^2 + eta*norm(as.matrix(mysidainner$myL[[d]])%*%myalphaold[[d]],'f')^2))
    relObj=abs(ObjNew-ObjOld)/ObjOld
    print(reldiff)
    print(diffalpha)
    print(relObj)
  }

  #classification
  if(AssignClassMethod=='Joint'){
  myclassify=sidaclassify(myalpha,Xtestdata,Xdata,Y,AssignClassMethod='Joint')
  sidaneterror=sum(myclassify$PredictedClass!=Ytest)/length(Ytest)
  }else if(AssignClassMethod=='Separate'){
  sidaneterror=sapply(1:(nK+1), function(x)  sum(myclassify$PredictedClass[,x]!=Ytest)/length(Ytest) )
  }


  #sum pairwise correlations of training data
  ss=list()
  #sum pairwise correlations
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xtestdata[[d]]%*%myalpha[[d]]
      X2=Xtestdata[[j]]%*%myalpha[[j]]
      X1=scale(X1, center=TRUE,scale=FALSE)
      X2=scale(X2, center=TRUE,scale=FALSE)
      X1X2=t(X1)%*%X2/dim(X1)[1]
      X1X1=t(X1)%*%X1/dim(X1)[1]
      X2X2=t(X2)%*%X2/dim(X2)[1]
      sumcorr3=sum(diag(X1X2%*%t(X1X2)))/(sqrt(sum(diag(X1X1%*%X1X1)))*sqrt(sum(diag(X2X2%*%X2X2))))
      sumCorr2=sumCorr2+sumcorr3
    }
    ss[[d]]=sumCorr2/length(dd)
  }

  sidanetcorrelation=sum(do.call(rbind,ss))/D


  result=list(sidaneterror=sidaneterror,sidanetcorrelation=sidanetcorrelation,hatalpha=myalpha,PredictedClass=myclassify$PredictedClass)
  return(result)

}
