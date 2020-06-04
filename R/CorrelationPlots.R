CorrelationPlots=function(Xtestdata,Ytest,hatalpha){

dsizes=lapply(Xtestdata, function(x) dim(x))
D=length(dsizes)
mycomb=t(combn(D, 2, FUN = NULL, simplify = TRUE))

nc=unique(Ytest)
nnc=length(unique(Ytest))
mycol=mat.or.vec(length(Ytest),1)
mypch=mat.or.vec(length(Ytest),1)
for(j in 1:length(mycol)){
  mycol[j]=nc[Ytest[j]]
  mypch[j]=nc[Ytest[j]]
}

for(d in 1:dim(mycomb)[1]){
  dd=mycomb[d,]
  Scoresd=Xtestdata[[dd[1]]]%*%hatalpha[[dd[1]]]
  Scoresj=Xtestdata[[dd[2]]]%*%hatalpha[[dd[2]]]
  Scores=cbind.data.frame(Scoresd[,1], Scoresj[,1], Ytest)
  colnames(Scores)=c("Disc1", "Disc2", "Class")

  #calculate RV coefficient
  X1=Xtestdata[[dd[1]]]%*%hatalpha[[dd[1]]]
  X2=Xtestdata[[dd[2]]]%*%hatalpha[[dd[2]]]
  X1=scale(X1, center=TRUE,scale=FALSE)
  X2=scale(X2, center=TRUE,scale=FALSE)
  X1X2=t(X1)%*%X2/dim(X1)[1]
  X1X1=t(X1)%*%X1/dim(X1)[1]
  X2X2=t(X2)%*%X2/dim(X2)[1]
  RVCoeff=sum(diag(X1X2%*%t(X1X2)))/(sqrt(sum(diag(X1X1%*%X1X1)))*sqrt(sum(diag(X2X2%*%X2X2))))
  RVCoeff=round(RVCoeff,digits=2)



  plot(Scores[,1], Scores[,2],col=mycol,lwd=3,pch=mypch
       ,xlab=paste(
         "First Discriminant Score for View", dd[1]),ylab=paste("Second Discriminant Score for View", dd[2]),xaxt="n",
       yaxt="n", main=paste("Correlation plot for views",dd[1], "and" ,dd[2], ",", "\u03C1 =", RVCoeff))
  legend("topleft",bty = "n",cex=0.8,legend=c(nc) ,col = 1:max(nc), pch=1:max(nc),title="Class")
}
return(NULL)


}
