cd("/Users/aliaghdaei/Desktop/Hypergraph/HypergraphFlowClustering-master/Exp-Amazon/Partitioning/")


using TickTock
using SparseArrays
using LinearAlgebra
using Clustering
using NearestNeighbors
using Distances
using Metis
using Laplacians
using Arpack
using Plots
using Statistics
using DelimitedFiles


include("hmet2ar.jl")
include("Star.jl")
include("Filter.jl")
include("h_score3.jl")
include("h_score.jl")
include("mx_func.jl")





ar = hmet2ar("ibm01.hgr")

mx = mx_func(ar)

A = Star(ar)

L = lap(A)

Nrv = 10

NumV = 5

SV = zeros(mx, Nrv * 5)

Qvec = zeros(Float64, 0)

iter = 0

for ii = 1:Nrv

    rv = (rand(Float64, size(A, 1), 5) .- 0.5).*2

    sm = zeros(mx, 5)

    for jj = 1:5

        sm[:, jj] = Filter(rv[:, jj], ii, A, mx)

        Q = h_score3(ar, sm[:, jj])

        append!(Qvec, sum(Q))

    end #end of jj

    SV[:, iter+1: iter + 5] = sm

    global iter += 5

end #end of for ii

Redge = zeros(Float64, 0)

for ii = 1:length(ar)

    nd = ar[ii]

    Rmx = 0

    for jj = 1:length(nd) - 1

        for kk = jj+1:length(nd)

            Rvec = zeros(Float64, 0)

            for ll = 1:size(SV,2)

                R = ((SV[nd[jj], ll] - SV[nd[kk], ll]) ^ 2)

                append!(Rvec, R)

            end #end of for ll

            Rvec = sort!(Rvec)

            Rnd = sum(Rvec[end - NumV + 1: end]) / NumV

            Rmx = max(Rnd, Rmx)

        end #end of for kk

    end #end of for jj

    append!(Redge, Rmx)

end #end of for ii

RedR = 0.1

N = round(Int, RedR * length(ar))

Pos = sortperm(Redge)

flag = falses(mx)

flagE = falses(length(ar))

val = 1

idx = zeros(Int, mx)

for ii = 1:N

    nd = ar[Pos[ii]]

    fg = flag[nd]

    fd1 = findall(x->x==1, fg)

    if length(fd1) == 0

        flagE[Pos[ii]] = 1

        idx[nd].= val

        flag[nd] .= 1

        global val +=1

    end # endof if

end #end of for ii

fdz = findall(x-> x==0, idx)

fdnz = findall(x-> x!=0, idx)

V = vec(val:val+length(fdz)-1)

idx[fdz] = V

fdE = findall(x->x==1, flagE)

deleteat!(ar, fdE)

ar_new = Any[]

for ii = 1:length(ar)

    nd = ar[ii]

    nd_new = unique(idx[nd])

    push!(ar_new, nd_new)

end #end of for ii



## compute the original non-linear quadratic value
ar = hmet2ar("ibm03.hgr")

A =Star(ar)

L = lap(A)

dg = sum(A, dims = 1) .^ (-.5)

I2 = 1:size(L,1)

D = sparse(I2, I2, sparsevec(dg))

L = D * L * D

L[diagind(L,0)] = L[diagind(L,0)] .+ 0.01

EV = eigs(L, nev = 3, which= :SM)

V = real(EV[2])

mx = mx_func(ar)

V = V[1:mx, 2]

Q_max = sum(h_score3(ar, V))

Q_mediator = sum(h_score(ar, V))

println("Original Q_max: ", Q_max)

println("Original Q_mediator: ", Q_mediator)

println("V_org: ", mx)

println("E_org: ", length(ar))

## compute the new matrix non-linear quadratic value

A =Star(ar_new)

L = lap(A)

dg = sum(A, dims = 1) .^ (-.5)

I2 = 1:size(L,1)

D = sparse(I2, I2, sparsevec(dg))

L = D * L * D

L[diagind(L,0)] = L[diagind(L,0)] .+ 0.01

EV = eigs(L, nev = 3, which= :SM)

V = real(EV[2])

mx_new = mx_func(ar_new)

V = V[1:mx_new, 2]

Q_max = sum(h_score3(ar_new, V))

Q_mediator = sum(h_score(ar_new, V))

println("New Q_max: ", Q_max)

println("New Q_mediator: ", Q_mediator)

println("V_new: ", mx_new)

println("E_new: ", length(ar_new))
