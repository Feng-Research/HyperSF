using SparseArrays
using MatrixNetworks
using LinearAlgebra
using Random

include("HyperSF.jl")
include("Functions.jl")
include("../include/HyperLocal.jl")
include("../include/Helper_Functions.jl")
include("../include/maxflow.jl")


filename = "ibm01.hgr"

## L controls the coarsening ratio by applying L-levels of k-mean clustering
L = 4

## R adjusts the ratio of selected clusters (low-quality clusters)
# for applying the flow-based technique (0<R<=1)
R = .1

## IDN is the index cluster that is assigned to each node
IDN = HyperSF(filename, L, R)



## Writing the coarse hypergraph into output file
cd("../data/")
ar = ReadInp(filename)
cd("../src/")
ar_new = Any[]
@inbounds for ii = 1:length(ar)
    
    nd = ar[ii]
    ndN = unique(IDN[nd])
    push!(ar_new, sort(ndN))

end #end of for ii

ar_new = unique(ar_new)

### removing hyperedges with cardinality of 1
HH = INC(ar_new)
ss = sum(HH, dims=2)
fd1 = findall(x->x==1, ss[:,1])
deleteat!(ar_new, fd1)


Whgr("Out.hgr", ar_new)
