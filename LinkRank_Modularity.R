# LinkRank (using power method)
lr.modularity <- function(g,partition, 
                          alpha=.85, 
                          max.iter=50000,
                          tol = .000001,
                          use.weights=T) {
                          
   ## Input :
   
   ## g           = graph (igraph object)
   ## partition   = graph partition (vector)
   ## alpha       = teleportation probability
   ## max.iter    = maximum iterations to calculate Perron vector
   ## tol         = tolerance
   ## use.weights = whether to use edge-weights in the calculation
   
   require(igraph)
   
   if (is.null(vertex_attr(g,'name'))) 
      stop('vertex attribute "name" is NULL')
   
   m.list <- cbind(vertex_attr(g,'name'),partition)
   
   if (use.weights) {
      if ('weight' %in% edge_attr_names(g)==F) {
         stop('No edge attribute named "weight" found for modularity calculation')
      }
      A <- get.adjacency(g, type='both', attr='weight',sparse=F)
   } else A <- get.adjacency(g, type='both', sparse=F)
   
   # no of nodes
   n <- nrow(A)
   # out degrees
   out.deg <- rowSums(A)
   # dead-end nodes
   dangling <- out.deg==0
   # row-normalize A
   G.temp <- sweep(A, 1, out.deg, FUN='/') 
   # set rows for dead-end nodes to zero
   G.temp[dangling,] <- 0
   # generate Google matrix
   Adj <- matrix(1/n * (alpha*dangling + 1- alpha), nrow=n, ncol=n)
   G <- alpha*G.temp + Adj
   # generate dominant eigenvalue and eigenvector by power method
   v.new <- runif(n)
   v.new <- v.new/sum(v.new)
   v.old <- rep(1,n)
   iter <- 0

   while (!isTRUE(all.equal(v.old,v.new,tolerance=tol)) & (iter <= max.iter)) {
      
      v.old <- v.new
      tmp.v <- v.old %*% G
      v.new <- tmp.v/sqrt(sum(tmp.v*tmp.v))
      e.val <- sum(v.old*v.new)/sum(v.new*v.new)
      
      if (iter==max.iter) {
         warning(paste0('Max.iter reached! Eigenvalue at iteration ',
                        max.iter, ' is ', e.val))
      }
      
      iter <- iter + 1
   }
   p.vec <- v.new/sum(v.new)
   
   # LinkRank matrix
   Q <- sweep(G,1,p.vec, '*') -  crossprod(p.vec)
   
   # LinkRank Modularity
   v.names <- vertex_attr(g,'name')
   parts <- sort(unique(partition))
   b.sum <- sapply(parts, function(z) {
      p.set <- v.names[partition==z]
      sum(Q[p.set,p.set])
   })
   
   return(sum(b.sum))
}
