R code to calculate LinkRank modularity (Kim et al. 2012). 

The LinkRank modularity is a generalization of modularity to directed and weighted networks. While Leicht and Newman (2008) has proposed an alternative modularity measure for weighted and directed networks, it has been pointed out that it does not appropriately incorporate the direction of the ties (see, for example, Fortunato 2010). The LinkRank modularity measure can be intuitively understood as the difference between the fraction of the time spend by a random walker within the communities and the expected value of that fraction under a null model in which the random walker "jumps" freely around given the constraint that the PageRank of each node is preserved. Kim et al. (2012) show that in the case of an undirected network, the LinkRank modularity reduces to the modularity proposed by Leicht and Newman (2008), so that the former can be considered as a generalization of the latter.

Notes:

1. The code uses the `page_rank` function from the `igraph` package to calculate the Perron vector. 
2. The default teleportation probability is set to .15 as recommended in Kim et al (2012). The option `damping` is equal to one minus the teleportation probability (in order to keep it inline with the `igraph::page_rank` function). Notice that the teleportation probability is often not needed in many undirected and connected graphs. Teleportation can be set off by specifying `damping = 1`.
3. Similar to the `weights` option in many `igraph` functions, if `weights` is NULL and the graph has a weight edge attribute then that is used. If weights is a numerical vector then it used, even if the graph has a weights edge attribute. If this is NA, then no edge weights are used (even if the graph has a weight edge attribute).
4. As the code is written in the `R` language, it is rather slow. So, the provided function might not suited for the detection of community structures in moderate to large networks. For the evaluation of community structures, it should run fine.

#### References

Fortunato. 2010. "Community Detection in Graphs," _Physics Reports_, 486

Kim, Y., S.W. Son, H. Jeong. 2010. "Finding Communities in Directed Networks," _Physical Review E_, 81

Leicht, E. A. and M. E. Newman. 2008. "Community Structure in Directed Netowrks," _Physical Review Letters_, 100
