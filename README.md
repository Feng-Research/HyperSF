# HyperSF
HyperSF: Spectral Hypergraph Coarsening via Flow-based Local Clustering (Accepted by ICCAD'21) 
Link to the paper: https://arxiv.org/abs/2108.07901)

Authors: Ali Aghdaei, Zhiqiang Zhao, Zhuo Feng

Abstract:
Hypergraphs allow modeling problems with multi-way high-order relationships. However, the computational cost of most existing hypergraph-based algorithms can be heavily dependent upon the input hypergraph sizes. To address the ever-increasing computational challenges, graph coarsening can be potentially applied for preprocessing a given hypergraph by aggressively aggregating its vertices (nodes). However, state-of-the-art hypergraph partitioning (clustering) methods that incorporate heuristic graph coarsening techniques are not optimized for preserving the structural (global) properties of hypergraphs. In this work, we propose an efficient spectral hypergraph coarsening scheme (HyperSF) for well preserving the original spectral (structural) properties of hypergraphs. Our approach leverages a recent strongly-local max-flow-based clustering algorithm for detecting the sets of hypergraph vertices that minimize ratio cut. To further improve the algorithm efficiency, we propose a divide-and-conquer scheme by leveraging spectral clustering of the bipartite graphs corresponding to the original hypergraphs. Our experimental results for a variety of hypergraphs extracted from real-world VLSI design benchmarks show that the proposed hypergraph coarsening algorithm can significantly improve the multi-way conductance of hypergraph clustering as well as runtime efficiency when compared with existing state-of-the-art algorithms.

Julia Version: 1.5.3

![Overview7 (1)](https://user-images.githubusercontent.com/85693952/128752511-64572dd2-aff6-4126-9290-c91a78b4c649.png)
