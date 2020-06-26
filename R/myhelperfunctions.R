myfastIDAnonsparse=function(Xdata, Y,weight){


  #   %--------------------------------------------------------------------------
  #   %myfastIDAnonsparse.R: function to obtain nonsparse solution to integrative lda problem
  # %and to obtain matrix needed in constraints
  # %--------------------------------------------------------------------------
  #
  #

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }


  #define weights
  w1=weight;
  w2=2*(1-weight)/(D*(D-1))

  Y=as.vector(Y)

  nc=max(unique(as.vector(Y)))

  Crxd=list()
  Sbx=list()
  myalphaold1=list()
  myalphaold2=list()
  myalphaoldmat=list()
  rmyalphaoldmat=list()
  sqrtminvmat=list()
  tildealphamat=list()
  tildelambda=list()

  for (d in 1:D){
    myX=as.matrix(Xdata[[d]])
    n=dim(myX)[1]
    p=dim(myX)[2]
    mysvd=svd(t(myX));
    Ux1=mysvd$u;
    V=mysvd$v;
    W=diag(mysvd$d)
    R=W%*%t(V)

    rdata2=cbind(Y,t(R))
    rdata=t(rdata2)
    mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
    mr=rowMeans(rdata[-1,])

    nc=max(unique(Y))
    C=list()
    for(i in 1:nc)
    {
      C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
    }
    C=as.matrix(do.call(rbind,C))
    Swrx=t(C)%*%C /(n-1)

    #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
    Crx=R-rowMeans(R)
    Srx=Crx%*%t(Crx)/(n-1)

    Sbrx=Srx-Swrx
    Sbrx=Sbrx + t(Sbrx)
    Sbx[[d]]=Sbrx

    lambda=sqrt(log(p)/n)
    if(n<p){
      Strx=Swrx + lambda*diag(n)
      Strx=Strx + t(Strx)

    } else if(n >= p){
      Strx=Swrx
      Strx=Strx + t(Strx)
    }

    Crxd[[d]]=Crx;

    #Set mybetaold and myalphaold as LDA solutions

    sqrtminv= mysqrtminv(Strx)$sqrtminv;
    sqrtminvmat[[d]]=sqrtminv;

    
    myeigen=eigs_sym(sqrtminv%*%Sbrx%*%sqrtminv,nc-1,which="LM")
    #myeigen=eigen(sqrtminv%*%Sbrx%*%sqrtminv,symmetric=TRUE)
    
    myalphaold1[[1]]=Ux1%*%myeigen$vectors
    myalphaoldmat[[d]]=do.call(rbind,lapply(myalphaold1, function(x) x/norm(x,'2')))
    rmyalphaoldmat[[d]]=myeigen$vectors
  }



  #nonsparse solution to integrative LDA
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #cross-covariance
    rSumassociation=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      myalphaold=rmyalphaoldmat[[j]];
      Sdj=Crxd[[d]]%*%t(Crxd[[j]])/(n-1);
      rassociation=Sdj%*%myalphaold%*%t(myalphaold)%*%t(Sdj)
      rSumassociation=rSumassociation + rassociation + t(rassociation);
    }
    #solution to integrative LDA
    myinteig=eigs_sym( sqrtminvmat[[d]]%*%( w1*Sbx[[d]] +  w2*rSumassociation)%*%sqrtminvmat[[d]],nc-1,which="LM");
    #myinteig=eigen(sqrtminvmat[[d]]%*%( w1*Sbx[[d]] +  w2*rSumassociation)%*%sqrtminvmat[[d]],symmetric=TRUE)
    myalphaold2[[1]]=myinteig$vectors
    tildealphamat[[d]]=do.call(rbind,lapply(myalphaold2, function(x) x/norm(x,'2')))
    tildelambda[[d]]=myinteig$values
  }

  result=list(tildealphamat=tildealphamat,tildelambda=tildelambda, myalphaoldmat=myalphaoldmat,sqrtminvmat=sqrtminvmat);
  return(result)


}

myfastinner = function(Xdata,Y,sqrtminv,myalphaoldmat,tildealphamat, weight=0.5){

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }


  #define weights
  w1=weight;
  w2=2*(1-weight)/(D*(D-1))

  Y=as.vector(Y)

  separationd=list()
  associationd=list()
  SepAndAssoc=list()
  SepAndAssocd=list()
  Crxd=list()
  Ux=list()

  for (d in 1:D){
    myX=as.matrix(Xdata[[d]])
    n=dim(myX)[1]
    p=dim(myX)[2]
    mysvd=svd(t(myX));
    Ux1=mysvd$u;
    V=mysvd$v;
    W=diag(mysvd$d)
    R=W%*%t(V)
    rdata2=cbind(Y,t(R))
    rdata=t(rdata2)
    mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
    mr=rowMeans(rdata[-1,])

    nc=length(unique(as.vector(Y)))
    C=list()
    for(i in 1:nc)
    {
      C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
    }
    C=as.matrix(do.call(rbind,C))
    Swrx=t(C)%*%C /(n-1)

    #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
    Crx=R-rowMeans(R)
    Srx=Crx%*%t(Crx)/(n-1)

    Sbrx=Srx-Swrx
    Sbrx=Sbrx + t(Sbrx)
    #Sbx[[d]]=Sbrx

    lambda=sqrt(log(p)/n)
    if(n<p){
      Strx=Swrx + lambda*diag(n)
      Strx=Strx + t(Strx)

    } else if(n >= p){
      Strx=Swrx
      Strx=Strx + t(Strx)
    }

    Ux[[d]]=Ux1;
    Crxd[[d]]=Crx;
    sqrtminvStrx = sqrtminv[[d]]
    
    separationd[[d]]=w1*Ux1%*%sqrtminvStrx%*%Sbrx%*%sqrtminvStrx%*%tildealphamat[[d]]
  }


  #obtain association matrices

  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #cross-covariance
    Sumassociation=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      myalphaold=myalphaoldmat[[j]]
      Sdj=Crxd[[d]]%*%t(Crxd[[j]])
      assoc2=(Sdj%*%t(Ux[[j]])%*%(myalphaold%*%t(myalphaold))%*%Ux[[j]]%*%t(Sdj));
      association=Ux[[d]]%*%sqrtminv[[d]]%*%(assoc2+t(assoc2))%*%(t(Ux[[d]])%*%Ux[[d]])%*%sqrtminv[[d]]%*%tildealphamat[[d]]/((n-1)^2);
      Sumassociation=Sumassociation + association ;
    }

    SepAndAssoc[[1]]=separationd[[d]]+ w2*Sumassociation;
    SepAndAssoc2=lapply(SepAndAssoc, function(x) x/norm(x,'i'));
    SepAndAssocd[[d]]= do.call(rbind,SepAndAssoc2)

  }

  result=list(SepAndAssocd=SepAndAssocd, Ux=Ux);
  return(result)
}

mysqrtminv=function(W){
  #W is symmetric, positive definite
  mysvd=svd(W);
  d=diag(mysvd$d^(-0.5))
  out=mysvd$u%*%d%*%t(mysvd$u)
  result=list(sqrtminv=out)
  return(result)
}

sidainner = function(Xdata,Y,sqrtminv,myalphaold,tildealpha, tildelambda,Tau,weight,withCov){


  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }

  #if withCov=TRUE, set tuning parameter as 0
  if(withCov==TRUE){
    Tau[[D]]=0.00001
  }

  nK=length(unique(as.vector(Y))) -1


  myfastnslda=myfastinner(Xdata,Y, sqrtminv, myalphaold,tildealpha, weight);

  myhatalpha=list()

  for(d in 1:D){
    p=dim(Xdata[[d]])[2]

    ##solve for SIDA directions
     Alphai=Variable(p,nK)
     Objx=sum(norm2(Alphai,axis=1))
    #
    # #defining the constraints
    constraints=list(norm_inf(sum_entries(abs(as.matrix(myfastnslda$SepAndAssocd[[d]]) - Alphai%*%as.matrix(diag(tildelambda[[d]],nrow=length(tildelambda[[d]])))), axis=1 ))<= Tau[[d]])
    prob=Problem(Minimize(Objx),constraints)
    result=solve(prob,solver="ECOS")
    alphai=result$getValue(Alphai)
    alphai[abs(alphai) <=10^-5]=0

    if(sum(sum(abs(alphai)))==0){
      hatalpha=alphai
    } else {
      hatalpha=alphai/matrix(rep(t(sqrt(colSums(alphai*alphai))),times=p),ncol=ncol(alphai),byrow=TRUE)
      Q1 = qr(hatalpha)
      hatalpha=qr.Q(Q1)
      hatalpha[abs(hatalpha) <=10^-5]=0
    }
    myhatalpha[[d]]=hatalpha;

  }
  result=list(hatalpha=myhatalpha)
  return(result)
}


sidanetinner = function(Xdata,Y,sqrtminv,myalphaold,tildealpha, tildelambda,Tau,weight,eta,myedges,myedgeweight,mynormLaplacianG=NULL,withCov){

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }


  if(withCov==TRUE){
    Tau[[D]]=0.00001
  }

  if(is.null(mynormLaplacianG)){
    mynormLaplacianG=myNLaplacianG(Xdata,myedges,myedgeweight)
  }
  
  nK=length(unique(as.vector(Y))) -1


  myfastnslda=myfastinner(Xdata,Y, sqrtminv, myalphaold,tildealpha, weight)
  

  myhatalpha=list()
 # myL=list()
  
  for(d in 1:D){
    p=dim(Xdata[[d]])[2]

    #if edge information is empty, then no group information
    #utilizes sida;
    ee=c(as.matrix((myedges[[d]]!=0)*1))
    if(max(ee)==1){
      L2=mynormLaplacianG[[d]]
      #solve for SIDANet directions
      Alphai=Variable(p,nK)
      LB=as.matrix(L2)%*%Alphai
      Objx=eta%*%sum(norm2(LB,axis=1)) +(1-eta)%*%sum(norm2(Alphai,axis=1))

      #defining the constraints
      constraints=list(norm_inf(sum_entries(abs(as.matrix(myfastnslda$SepAndAssocd[[d]]) - Alphai%*%as.matrix(diag(tildelambda[[d]],nrow=length(tildelambda[[d]])))), axis=1 ))<= Tau[[d]])
      prob=Problem(Minimize(Objx),constraints)
      result=solve(prob,solver="ECOS")
      alphai=result$getValue(Alphai)
      alphai[abs(alphai) <=10^-5]=0
    }else if(max(myedges[[d]])==0){
       #if edge information is not available, use SIDA
      #solve for SIDANet directions
      Alphai=Variable(p,nK)
      Objx=sum(norm2(Alphai,axis=1))

      #defining the constraints
      constraints=list(norm_inf(sum_entries(abs(as.matrix(myfastnslda$SepAndAssocd[[d]]) - Alphai%*%as.matrix(diag(tildelambda[[d]],nrow=length(tildelambda[[d]])))), axis=1 ))<= Tau[[d]])
      prob=Problem(Minimize(Objx),constraints)
      result=solve(prob,solver="ECOS")
      alphai=result$getValue(Alphai)
      alphai[abs(alphai) <=10^-5]=0

    }


    #normalize alphai
    if(sum(sum(abs(alphai)))==0){
      hatalpha=alphai
    } else {
      hatalpha=alphai/matrix(rep(t(sqrt(colSums(alphai*alphai))),times=p),ncol=ncol(alphai),byrow=TRUE)
      Q1 = qr(hatalpha)
      hatalpha=qr.Q(Q1)
      hatalpha[abs(hatalpha) <=10^-5]=0
    }
    myhatalpha[[d]]=hatalpha;

  }
  result=list(hatalpha=myhatalpha,myL=mynormLaplacianG)
  return(result)
}

myNLaplacianG=function(Xdata,myedges,myedgeweight){
  myL=list()
  D = length(Xdata)
  for(d in 1:D){
    p=dim(Xdata[[d]])[2] 
    ee=c(as.matrix((myedges[[d]]!=0)*1))
    if(max(ee)==1){
      edgesd=as.matrix(myedges[[d]])
      edgeNode=unique(c(as.matrix(edgesd)))
      edgeweightd=myedgeweight[[d]]
      
      #create undirected graph from edge list
      G=graph_from_edgelist(edgesd, directed=FALSE)
      L2=Matrix(0, nrow=p,ncol=p,sparse=TRUE)
      
      #calculate normalized laplacian of the weighted or unweighted graph
      if(max(edgeweightd)!=0){
        nL=laplacian_matrix(G, normalized=TRUE, weights=edgeweightd, sparse=igraph_opt("sparsematrices"))
      }else if(max(edgeweightd)==0){
        nL=laplacian_matrix(G, normalized=TRUE, weights=NULL,sparse=igraph_opt("sparsematrices"))
      }
      L2[edgeNode,edgeNode]=nL
      myL[[d]]=Matrix(L2,sparse = TRUE)
    }else if(max(myedges[[d]])==0){
      myL[[d]]=Matrix(diag(rep(1,p)),sparse = TRUE)
    }
  }
  return(myL)
}

