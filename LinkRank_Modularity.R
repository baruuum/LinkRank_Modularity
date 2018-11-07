# LinkRank 
lr.modularity <- function(g,
                          partition, 
                          damping = .85, 
                          pr.algo = 'prpack',
                          weights = NULL) {
    
    ## g           = graph (igraph object)
    ## partition   = graph partition (numeric vector of memberships or "communities" object)
    ## damping     = damping factor (1 - teleportation prob.)
    ## pr.algo     = algorithm to calculate Perron vector,
    ##               possible options are "prpack", "arpack", and "power"
    ## weights     = If this is NULL and the graph has a weight edge attribute
    ##               then that is used. If weights is a numerical vector then
    ##               it used, even if the graph has a weights edge attribute.
    ##               If this is NA, then no edge weights are used (even if the
    ##               graph has a weight edge attribute)
    
    # check args
    if (!is.igraph(g)) 
        stop('graph is not an i.graph object')
    
    if (damping > 1 | damping < 0) 
        stop('damping factor has to be between zero and one!')
    
    # get algorithm name to calculate Perron vector
    pr.algo <- match.arg(pr.algo, c('prpack','arpack','power'))
    
    # no of nodes
    n <- vcount(g)
    # node sequence
    v.seq <- seq_len(n)
    
    # get membership vector
    if (class(partition) == 'communities') {
        
        pp <- membership(partition)
        
    } else {
      
      if (!is.numeric(partition))
          stop("'partition' has to be a 'communities' object or a numeric vector!")
      pp <- partition
      
    }
    
    # check dimensions
    if (length(pp) != n) 
        stop('Length of membership vector differs from number of nodes!')
    
    # get adjacency matrix & out-degree
    if (is.vector(weights) & length(weights) > 1) {
        
        # check args
        if (ecount(g) != length(weights))
            stop("'weights' differes in length from ecount!")
        if (!is.numeric(weights))
            stop("'weights' must be 'NA','NULL', or a numeric vector!")
        
        edge_attr(g, 'tmp') <- weights
        A <- get.adjacency(g, type = 'both', attr = 'tmp')
        
        out.deg <- strength(g, mode = 'out', weights = weights)
        
    } else if (is.null(weights)) {
        
        if ('weight' %in% edge_attr_names(g)) {
            
            A <- get.adjacency(g, type='both', attr='weight')
            out.deg <- strength(g, mode = 'out')
            
        }  else {
            
            A <- get.adjacency(g, type='both')
            out.deg <- degree(g, mode = 'out')
            
        }
        
    } else if (is.na(weights)) {
        
        A <- get.adjacency(g, type='both')
        out.deg <- degree(g, mode = 'out')
        
    } else {
        
        stop("'weights' option has to be 'NA','NULL', or a numeric vector!")
        
    }
    
    # dead-end nodes
    dangling <- out.deg == 0
    
    # row-normalize A (recycle vector)
    G.temp <- A / out.deg
    # equivalent to sweep(A, 1, out.deg, FUN='/')
    
    # set rows for dead-end nodes to zero
    if (sum(dangling) > 0) {
        G.temp[dangling,] <- 0
    }
    
    # add teleportation probabilities
    Tmat <- Matrix::Matrix(1/n * (damping * dangling + 1 - damping), 
                           nrow = n, ncol = n)
    G <- damping * G.temp + Tmat
    
    # get Perron vector (PageRank)
    p.vec <- page_rank(g, damping = damping, algo = pr.algo, weights = weights)$vector
    
    # LinkRank matrix
    Q <- G * p.vec -  tcrossprod(p.vec)
    # equivalent to sweep(G, 1, p.vec, '*') -  tcrossprod(p.vec)
    
    # get LinkRank Modularity by summing over within-community weights
    return(sum(Q[outer(pp, pp, '==')]))
  
}
