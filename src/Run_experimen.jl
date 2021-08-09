using SparseArrays
using LinearAlgebra
using Clustering
using Distances
using Metis
using Laplacians


include("HyperSF.jl")

Input = "../data/ibm01.hgr"

#Threshold on the conductance of the selected clusters; the lower, the better
#should be less than 1
CndT = 0.8

# AC is the average conductance of all clusters
AC, ar_coarse = HyperSF(Input, CndT)

# The incidence matrix corresponding to the coarsened hypergraph
H = INC3(ar_coarse)
