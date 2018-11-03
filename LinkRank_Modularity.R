# LinkRank 
lr.modularity <- function(g,
                          partition, 
                          damping = .85, 
                          pr.algo = 'prpack',
                          use.weights = T) {
    
    ## Input :
    
    ## g           = graph (igraph object)
    ## partition   = graph partition (vector or "communities" object)
    ## damping     = damping factor (1 - teleportation prob.)
    ## pr.algo     = algorithm to calculate Perron vector,
    ##               possible options are "prpack", "arpack", and "power"
    ## use.weights = whether to use edge-weights in the calculation
    
    require(igraph)
    
    # no of nodes
    n <- vcount(g)
    
    # node sequence
    v.seq <- seq_len(n)
    
    # get membership vector
    if (class(partition) == 'communities') {
        
        pp <- membership(partition)
        
    } else {
        
        pp <- partition
        
    } 
    
    # check dimensions
    if (length(pp) != n) 
        stop('Length of membership vector differs from number of nodes!')
    
    # get adjacency matrix
    if (use.weights) {
        
        if ('weight' %in% edge_attr_names(g) == F) {
            
            stop('No edge attribute named "weight" found!')
        }
        
        A <- get.adjacency(g, type='both', attr='weight')
        ww <- edge_attr(g, 'weight')
        
    } else {
        
        A <- get.adjacency(g, type='both')
        ww <- NULL
        
    }
    
    # out degrees
    if (use.weights) {
        
        out.deg <- strength(g, mode = 'out', weights = ww)
        
    } else {
        
        out.deg <- degree(g, mode = 'out')
        
    }
    
    # dead-end nodes
    dangling <- out.deg == 0
    
    # row-normalize A
    G.temp <- sweep(A, 1, out.deg, FUN='/')
    
    # set rows for dead-end nodes to zero
    if (sum(dangling) > 0) {
        
        G.temp[dangling,] <- 0
        
    }
    
    # add teleportation probabilities
    Tmat <- Matrix::Matrix(1/n * (damping * dangling + 1 - damping), 
                           nrow = n, ncol = n)
    G <- damping * G.temp + Tmat
    
    # get Perron vector (PageRank)
    p.vec <- page_rank(g, damping = damping, algo = pr.algo, weights = ww)$vector
    
    # LinkRank matrix
    Q <- sweep(G, 1, p.vec, '*') -  tcrossprod(p.vec)
    
    # get LinkRank Modularity
    parts <- sort(unique(pp))
    b.sum <- sapply(parts, function(z) {
        p.set <- v.seq[pp == z]
        sum(Q[p.set, p.set])
    })
    
    return(sum(b.sum))
}
