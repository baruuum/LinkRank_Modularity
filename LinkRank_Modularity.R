# LinkRank 
lr.modularity <- function(g,
                          partition, 
                          alpha=.85, 
                          pr.algo = 'prpack',
                          use.weights=T) {
    
    ## Input :
    
    ## g           = graph (igraph object)
    ## partition   = graph partition (vector or "communities" object)
    ## alpha       = damping factor (1 - teleportation prob.)
    ## pr.algo     = algorithm to calculate Perron vector 
    ## use.weights = whether to use edge-weights in the calculation
    
    require(igraph)
    
    # get vertex sequences
    
    if (is.null(vertex_attr(g, 'name'))) {
        
        warning('vertex attribute "name" is NULL ... using numerical values')
        v.names <- 1:vcount(g)
        
    } else {
        
        v.names <- vertex_attr(g, 'name')
        
    }
    
    # get membership vector
    
    if (class(partition) == 'communities') {
        
        pp <- membership(pp)
        
    } else {
        
        pp <- partition
        
    } 
    
    # check dimensions
    if (length(pp) != length(v.names)) 
        stop('Length of partition differs from number of nodes!')
    
    # matrix of vertices and partition
    m.mat <- cbind(v.names,partition)
    
    # get adjacency matrix
    if (use.weights) {
        
        if ('weight' %in% edge_attr_names(g)==F) {
            
            stop('No edge attribute named "weight" found!')
        }
        
        A <- get.adjacency(g, type='both', attr='weight')
        ww <- edge_attr(g, 'weight')
        
    } else {
        
        A <- get.adjacency(g, type='both')
        ww <- NULL
        
    }
    
    # no of nodes
    n <- vcount(g)
    
    # out degrees
    if (use.weights) {
        
        out.deg <- strength(g, mode = 'out', weights = ww)
        
    } else {
        
        out.deg <- degree(g, mode = 'out')
        
    }
    
    # dead-end nodes
    dangling <- out.deg==0
    
    # row-normalize A
    G.temp <- sweep(A, 1, out.deg, FUN='/')
    
    # set rows for dead-end nodes to zero
    if (sum(dangling) > 0) {
        
        G.temp[dangling,] <- 0
        
    }
    
    # add teleportation probabilities
    Tmat <- Matrix::Matrix(1/n * (alpha*dangling + 1- alpha), nrow=n, ncol=n)
    G <- alpha*G.temp + Tmat
    
    # get Perron vector (PageRank)
    p.vec <- page_rank(g, damping = alpha, algo = pr.algo, weights = ww)$vector
    
    # LinkRank matrix
    Q <- sweep(G,1,p.vec, '*') -  tcrossprod(p.vec)
    
    # get LinkRank Modularity
    parts <- sort(unique(partition))
    b.sum <- sapply(parts, function(z) {
        p.set <- v.names[partition==z]
        sum(Q[p.set,p.set])
    })
    
    return(sum(b.sum))
}
