sidaclassify=function(hatalpha,Xtestdata,Xdata,Y,AssignClassMethod='Joint', standardize=TRUE){

#hatalpha is a list of d estimated SIDA vectors, for each view
#Y is a vector of training observations

  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }

  D=length(Xtestdata)
  #classification
  Projtest=sapply(1:D, function(i) list(Xtestdata[[i]]%*%hatalpha[[i]]))
  Projtrain=sapply(1:D, function(i) list(Xdata[[i]]%*%hatalpha[[i]]))


nc= length(unique(as.vector(Y)))
ntest=dim(Projtest[[1]])[1]

PredclassSeparate=list()
if(AssignClassMethod=='Separate'){
  for(d in 1:D){
    ProjXtestdatad=Projtest[[d]]
    ProjXdatad=Projtrain[[d]]
    ProjXdata1d=cbind(Y,ProjXdatad)
    Projmv=aggregate(ProjXdata1d[,-1],list(ProjXdata1d[,1]),mean)
    distv=list()
    jrep=list()
    for(j in 1: nc){
      rProjm=matrix( rep(Projmv[j,-1],times= ntest), ncol=ncol(ProjXdatad), byrow=TRUE)
      #euclidean distance
      sqdiff=(ProjXtestdatad-as.numeric(rProjm))^2
      dist1=rowSums(sqdiff)^0.5
      jrep[[j]]=j*rep(1,times=ntest)
      distv[[j]]=dist1
    }
    distv=do.call(cbind,distv)
    dim(distv)=c(nc*nrow(distv),1)

    jrep=do.call(cbind,jrep)
    dim(jrep)=c(nc*nrow(jrep),1)

    distv=cbind(jrep, distv)

    #The following code outputs the assigned class
    rdistvX1=matrix(distv[,-1], nrow=ntest,ncol=nc)
    minX1=apply(rdistvX1,1,min)
    minind=which(rdistvX1==minX1, arr.ind=TRUE) #minimum indices
    predclassX1=minind[order(minind[,1]),2]
    PredclassSeparate[[d]]=predclassX1
  }
  Predclass=do.call(cbind,PredclassSeparate)
}
else if(AssignClassMethod=='Joint') {
  #classification for joint
  ProjtestJoint=do.call(cbind,Projtest)
  ProjtrainJoint=do.call(cbind, Projtrain)

  ntest=dim(ProjtestJoint)[1]
  ProjtrainJointX1=cbind(Y,ProjtrainJoint)
  Projmv=aggregate(ProjtrainJointX1[,-1],list(ProjtrainJointX1[,1]),mean)

  distv=list()
  jrep=list()
  for(j in 1: nc){
    rProjm=matrix( rep(Projmv[j,-1],times= ntest), ncol=ncol(ProjtrainJoint), byrow=TRUE)
    #euclidean distance
    sqdiff=(ProjtestJoint-as.numeric(rProjm))^2
    dist1=rowSums(sqdiff)^0.5
    jrep[[j]]=j*rep(1,times=ntest)
    distv[[j]]=dist1
  }
  distv=do.call(cbind,distv)
  dim(distv)=c(nc*nrow(distv),1)

  jrep=do.call(cbind,jrep)
  dim(jrep)=c(nc*nrow(jrep),1)

  distv=cbind(jrep, distv)


  #The following code outputs the assigned class
  rdistvX1=as.data.frame(matrix(distv[,-1], nrow=ntest,ncol=nc))
  minX1=apply(as.matrix(rdistvX1),1,min)
  minind=which(rdistvX1==minX1, arr.ind=TRUE) #minimum indices
  Predclass=minind[order(minind[,1]),2]
#  Predclass=predclassJ
}
result=list(PredictedClass=Predclass, AssignClassMethod=AssignClassMethod)
return(result)
}
