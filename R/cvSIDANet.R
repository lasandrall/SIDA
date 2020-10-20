cvSIDANet=function(Xdata=Xdata,Y=Y,myedges=myedges,myedgeweight=myedgeweight,withCov=FALSE,plotIt=FALSE,Xtestdata=NULL,Ytest=NULL,isParallel=TRUE,ncores=NULL,gridMethod='RandomSearch',AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,maxiteration=20, weight=0.5,thresh=1e-03,eta=0.5){


  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
    if(D==1){
      stop("There should be at least two datasets")
    }
  } else {
    stop("Input data should be a list")
  }


  #set defaults
  #If testing data are not provided, the default is to use training data
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
  }

  if(is.null(Ytest)){
    Ytest=Y
  }

  #check inputs for training data
  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  nsizes=lapply(Xdata, function(x) dim(x)[1])

  if(all(nsizes!=nsizes[[1]])){
    stop('The  datasets  have different number of observations')
  }

  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])
  if(all(ntestsizes!=ntestsizes[[1]])){
    stop('The testing datasets  have different number of observations')
  }


  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
    if(D==1){
      stop("There should be at least two datasets")
    }
  } else {
    stop("Input data should be a list")
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
  
  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }
  
  if(is.null(gridMethod)){
    gridMethod='RandomSearch'
  }
  
  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }
  
  if(is.null(isParallel)){
    isParallel=TRUE
  }
  
  if(is.null(nfolds)){
    nfolds=5
  }
  
  if(is.null(ngrid)){
    ngrid=8
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
  
  if(is.null(eta)){
    eta=1e-03
  }
  

  set.seed(1234)
  nK=length(unique(as.vector(Y))) -1

  nc=length(unique(as.vector(Y)))
  Nn=mat.or.vec(nc,1)
  foldid=list()
  for(i in 1:nc)
  {
    Nn[i]=sum(Y==i)
    mod1=Nn[i]%%nfolds
    if(mod1==0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds)),Nn[i])
    }else if(mod1> 0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds), 1:(Nn[i]%%nfolds)),Nn[i])
    }
  }

  foldid=unlist(foldid)

  #obtain tuning range common to all K


  if(withCov==TRUE){
    Dnew=D-1
  }else if(withCov==FALSE){
    Dnew=D
  }

  if(Dnew>2){
    ngrid=5
  }

  #calculate the normalized laplacian
  mynormLaplacianG=myNLaplacianG(Xdata=Xdata,myedges=myedges,myedgeweight=myedgeweight)
  myTauvec=sidanettunerange(Xdata,Y,ngrid,standardize,weight,eta,myedges,myedgeweight,withCov)

  #define the grid
  mygrid=expand.grid(do.call(cbind,myTauvec))
  gridcomb=dim(mygrid)[1]
  if(gridMethod=='RandomSearch'){
    if(Dnew<2){
      ntrials=ngrid}
      else if(Dnew==2){
      ntrials=floor(0.2*gridcomb)}
      else if(Dnew>2){
        ntrials=floor(0.15*gridcomb)
      }
    mytune=sample(1:gridcomb, ntrials, replace = FALSE)
    gridValues=mygrid[mytune,]
    }else if(gridMethod=='GridSearch'){
    gridValues=mygrid
  }


#  start_time=Sys.time()
  gridValues=as.matrix(gridValues)
  CVOut=matrix(0, nfolds, nrow(as.matrix(gridValues)))
  #cross validation
if(isParallel==TRUE){
  registerDoParallel()
  if(is.null(ncores)){
    ncores=parallel::detectCores()
    ncores=ceiling(ncores/2)}
  cl=makeCluster(ncores)
  #start_time=Sys.time()
  registerDoParallel(cl)
  CVOut=matrix(0, nrow(gridValues), nfolds)
  mycv=foreach(i = 1:nrow(gridValues), .combine='rbind',.export=c('sidanet','sidanetinner','myfastinner','myfastIDAnonsparse','mysqrtminv','sidaclassify', 'sidanettunerange','DiscriminantPlots','CorrelationPlots'),.packages=c('CVXR','RSpectra','igraph','Matrix')) %dopar% {
    mTau=sapply(1:D, function(itau) list(t(gridValues[,itau][i])))
    start_time=Sys.time()
    CVOut[i,]= sapply(1:nfolds, function(j){
      testInd=which(foldid==j)
      testX=lapply(Xdata, function(x) x[testInd,])
      testY=Y[testInd]
      trainX=lapply(Xdata, function(x) x[-testInd,])
      trainY=Y[-testInd]
      mysida=sidanet(Xdata=trainX,Y=trainY,myedges,myedgeweight,Tau=mTau,withCov,Xtestdata=testX,Ytest=testY,AssignClassMethod,plotIt=NULL,standardize,maxiteration,weight,thresh,eta,mynormLaplacianG)
      return(min(mysida$sidaneterror))
    } )
  }
  #end_time=Sys.time()      
  #SIDAtime=end_time-start_time
  CVOut=t(mycv)
  stopCluster(cl)
}else if(isParallel==FALSE){
  CVOut=matrix(0, nfolds, nrow(gridValues))
  for (i in 1:nfolds){
    testInd=which(foldid==i)
    testX=lapply(Xdata, function(x) x[testInd,])
    testY=Y[testInd]
    trainX=lapply(Xdata, function(x) x[-testInd,])
    trainY=Y[-testInd]

    cat("Begin CV-fold", i, "\n")
    CVOut[i,]= sapply(1:nrow(gridValues), function(itau){
    mTau=sapply(1:D, function(d) list(t(gridValues[itau,][d])))
    mysida=sidanet(Xdata=trainX,Y=trainY,myedges,myedgeweight,Tau=mTau,withCov,Xtestdata=testX,Ytest=testY,AssignClassMethod,plotIt=NULL,standardize,maxiteration,weight,thresh,eta,mynormLaplacianG)
      return(min(mysida$sidaneterror))
     } )
  }
}


  #compute average classification error
  minEorrInd=max(which(colMeans(CVOut)==min(colMeans(CVOut))))
  optTau=gridValues[ minEorrInd,]

  #Apply on testing data
  moptTau=sapply(1:D, function(i) list(t(gridValues[minEorrInd,][i])))
  mysida=sidanet(Xdata=Xdata,Y=Y,myedges,myedgeweight,Tau=moptTau,withCov,Xtestdata=Xtestdata,Ytest=Ytest,AssignClassMethod,plotIt=NULL,standardize,maxiteration,weight,thresh,eta,mynormLaplacianG)


  ss=list()
  #sum pairwise RV coefficients
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xtestdata[[d]]%*%mysida$hatalpha[[d]]
      X2=Xtestdata[[j]]%*%mysida$hatalpha[[j]]
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

  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    DiscriminantPlots(Xtestdata,Ytest,mysida$hatalpha)
    CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)
  }else{
    myDiscPlot=NULL
    myCorrPlot=NULL
  }

  result=list(CVOut=CVOut,sidaneterror=mysida$sidaneterror,sidanetcorrelation=sidanetcorrelation,hatalpha=mysida$hatalpha,PredictedClass=mysida$PredictedClass, optTau=moptTau,gridValues=gridValues, AssignClassMethod=AssignClassMethod, gridMethod=gridMethod)
  return(result)
}
