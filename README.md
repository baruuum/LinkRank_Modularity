R code to calculate LinkRank modularity (Kim et al. 2012). 

The LinkRank modularity is a generalization of modularity to directed and weighted networks. While Leicht and Newman (2008) has proposed an alternative modularity measure for weighted and directed networks, it has been pointed out that it does not appropriately incorporate the direction of the ties (see, for example, Fortunato 2010). The LinkRank modularity measure can be intuitively understood as the difference between the fraction of the time spend by a random walker within the communities and the expected value of that fraction under a null model in which the random walker "jumps" freely around given the constraint that the PageRank of each node is preserved. Kim et al. (2012) show that in the case of an undirected network, the LinkRank modularity reduces to the "directed" modularity proposed by Leicht and Newman (2008), so that the former can be considered as a generalization of the latter.

Note: This is not the most efficient code and probably not well-suited to analyze large graphs.

#### References

Fortunato. 2010. "Community Detection in Graphs," _Physics Reports_, 486

Kim, Y., S.W. Son, H. Jeong. 2010. "Finding Communities in Directed Networks," _Physical Review E_, 81

Leicht, E. A. and M. E. Newman. 2008. "Community Structure in Directed Netowrks," _Physical Review Letters_, 100
