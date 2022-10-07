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
       myL[[d]]=Matrix(diag(rep(1,p)))
     }
  }
  return(myL)
}