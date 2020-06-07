DiscriminantPlots=function(Xtestdata=Xtestdata,Ytest=Ytest,hatalpha=hatalpha){


  dsizes=lapply(Xtestdata, function(x) dim(x))
  D=length(dsizes)

  nc=unique(Ytest)
  if(length(unique(Ytest))==2){
    par(mfrow=c(1,D))
    for(d in 1:D){
      myScores=Xtestdata[[d]]%*%hatalpha[[d]]
      ss21=myScores[Ytest==nc[1],]
      ss22=myScores[Ytest==nc[2],]
      dd1=density(ss21);
      dd2=density(ss22);
      plot(dd1,xlim=c(min(myScores[,1]-0.7),max(myScores[,1]+0.7)),cex.axis=1.5,
           cex.lab=1.5,lwd=2.5,xlab="First Discriminant Score ",main=paste("Discriminant Plot for View",d),ylim=c(0,max(dd1$y,dd2$y)))
      lines(dd2,xlim=c(min(myScores[,1]-0.7),max(myScores[,1]+0.7)),col="red",cex.axis=1.5,cex.lab=1.5,lwd=2.5)
      points(cbind(ss21,rep(0,length(ss21))),col="black",pch=3)
      points(cbind(ss22,rep(0,length(ss22))),col="red",pch=1)
      #legend("topleft", legend = c(paste("Class", nc[1]),paste("Class", nc[2])), col = c("red","black"),cex=1.2,lty=1,lwd=2.5)
      legend("topleft", inset=c(0,0),legend=c(nc), col = c("black","red"),pch=c(3,1),title="Class")
      
    }
    par(mfrow=c(1,1))
  }else if(length(unique(Ytest))>2){
    par(mfrow=c(1,D))
    nnc=length(unique(Ytest))
    mycol=mat.or.vec(length(Ytest),1)
    mypch=mat.or.vec(length(Ytest),1)
    for(j in 1:length(mycol)){
      mycol[j]=nc[Ytest[j]]
      mypch[j]=nc[Ytest[j]]
      
    }
    for(d in 1:D){
      Scores=cbind.data.frame(Ytest,Xtestdata[[d]]%*%cbind(hatalpha[[d]],hatalpha[[d]]))
      plot(Scores[,2], Scores[,3],col=mycol,lwd=2.5,pch=mypch
           ,xlab=paste(
             "First Discriminant Score for View", d),ylab=paste("Second Discriminant Score for View", d),xaxt="n",yaxt="n", main=paste("Discriminant Plot for View",d))
      par(xpd=TRUE)
      legend("topright",bty = "n",legend=c(nc),col = 1:max(nc), pch=1:max(nc),title="Class")
    }
    par(mfrow=c(1,1))
  }


  return(NULL)

}

