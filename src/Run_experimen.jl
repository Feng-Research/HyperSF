using SparseArrays
using LinearAlgebra
using Clustering
using Distances
using Metis
using Laplacians


include("HyperSF.jl")

Input = "../data/ibm01.hgr"

#Threshold on the conductance of the selected clusters; the lower, the better
CndT = 0.3

# AC is the average conductance of all clusters
AC = HyperSF(Input, CndT)
