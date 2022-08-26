# HyperSF
Spectral Hypergraph Coarsening Code for

HyperSF: Spectral Hypergraph Coarsening via Flow-based Local Clustering 
https://ieeexplore.ieee.org/document/9643555/metrics#metrics

Ali Aghdaei, Zhiqiang Zhao, Zhuo Feng

This is the improved version of HyperSF incorporating effective resistance clustering method.

![Overview7 (1)](https://user-images.githubusercontent.com/85693952/128752511-64572dd2-aff6-4126-9290-c91a78b4c649.png)

# Requirements
Julia Version: 1.5.3

Packages:

SparseArrays

LinearAlgebra

MatrixNetworks v1.0.1

RandomV06 v0.0.2

# Testing 
Run "Run_experiment.jl" to generate the coarsened hypergraph for the availble test cases.

# Input
The input format is in the format of hMetis that is every line corresponds to each hyperedge.
The first line: #hyperedges, #nodes.

L: is an integer (L>0) to adjust the initial coarsening ratio before applying the flow-based technique 

R: the ratio of selected clusters (low-quality clusters) for applying the flow-based technique (0<R<=1).

# Output
IDN: returns the assigned cluster to each node.

The code generate the coarsened hypergraph as "Out.hgr" as the same format of input hypergraph.
