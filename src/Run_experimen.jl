cd("/Users/aliaghdaei/Desktop/Hypergraph/HypergraphFlowClustering-master/Exp-Amazon/Partitioning/HyperSF/Github/")

using SparseArrays
using LinearAlgebra
using Clustering
using Distances
using Metis
using Laplacians


include("HyperSF.jl")
include("Functions.jl")
include("HyperLocal.jl")

Input = "ibm01.hgr"

#Threshold on the conductance of the selected clusters; the lower, the better
CndT = 0.3

# AC is the average conductance of all clusters
AC = HyperSF(Input, CndT)
