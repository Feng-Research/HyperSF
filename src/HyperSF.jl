include("Functions.jl")
include("../include/HyperLocal.jl")
include("../include/Helper_Functions.jl")
include("../include/maxflow.jl")

function HyperSF(Inp::Strings, CndT::Float64)


    avg_cond = zeros(Float64, 0)
    kway_cond = zeros(Float64, 0)
    sz_mat = zeros(Int, 0)
    RR = zeros(Float64, 0)
    ## parameters
    grow_lvl = 2

    grownum = 100
    k = 10
    metis_P = 100


    Cvec_lvl = Any[]

    idx_mat = Any[]

    ar = ReadInp(Inp)

    NH = HyperNodes(ar)

    H = INC3(ar)

    mx = mx_func(ar)

    MM = Star(ar)


    ## generate 10 random vectors
    RD = (rand(Float64, size(MM,1), 10) .- 0.5).*2

    SV = Filter(RD, k, MM, mx)

    Hscore = h_score3(ar ,SV)

    idx = zeros(Int, mx)

    idx_new = he_cluster14(idx, ar, Hscore, SV, mx)

    ## Applying Metis
    N = round(Int, mx/metis_P)

    P = Metis.partition(MM, N)
    P = P[1:mx]

    flag = falses(mx)

    idx_coarse = zeros(Int, mx)

    val = 0

    Cvec_FB = zeros(Float64, 0)

    cut_mat = zeros(Int, 0)

    vol_mat = zeros(Int, 0)


    Threads.@threads for loop = 1:N

        nd_met = findall(x->x==loop, P)

        ## finding the clusters each node belongs to
        CL = idx_new[nd_met]

        ## finding non-unique cluster
        NU = nonunique2!(CL)

        for ii = 1:length(NU)

            fd1 = findall(x->x==NU[ii], CL)

            nd1 = nd_met[fd1]

            seedN = nd1[.!flag[nd1]]

            if length(seedN)>1

                ## expanding the network around the seed nodes
                nd = copy(seedN)

                HE = Any[]

                for ee = 1:grow_lvl

                    HE = Any[]

                    for dd = 1:length(nd)

                        append!(HE, NH[nd[dd]])

                    end

                    HE = sort(unique(HE))

                    new_nd = Any[]

                    for mm=1:length(HE)

                        nnd = ar[HE[mm]]

                        append!(new_nd, nnd)

                    end #end of mm

                    nd = sort(unique(new_nd))

                end #end of grow_lvl

                seedN2 = findall(x->in(x, seedN), nd)


                IM = H[HE, nd]
                #println("size(IM): ", size(IM))

                IMt = sparse(IM')

                ## d is node degree
                d = vec(sum(IM,dims=1))

                epsilon = 1.0

                delta = 1.0

                order = round.(Int, sum(IM, dims = 2))
                order = order[:,1]


                OneHop = get_immediate_neighbors(IM,IMt,seedN2)
                Rmore = BestNeighbors(IM,d,seedN2,OneHop,grownum)
                R = union(Rmore,seedN2)
                #Rs = findall(x->in(x,seedN),R)  #Force
                Rs = findall(x->in(x,seedN2),R)  #Force

                S, lcond = HyperLocal(IM,IMt,order,d,R,epsilon,delta,Rs,true)

                volA = sum(d)

                S_org = nd[S]

                S_met = findall(x->in(x, S_org), nd_met)

                S_org = nd_met[S_met]

                ## finding the not discovered nodes
                fgs = flag[S_org]

                S_org2 = S_org[.!fgs]

                if length(S_org2) > 0

                    S3 = findall(x->in(x, S_org2), nd)

                    condR,volR, cutR = tl_cond(IM,S3,d,delta,volA,order)

                    if condR < CndT

                        append!(Cvec_FB, condR)

                        idx_coarse[S_org2] .= val

                        flag[S_org2] .= 1

                        val+=1

                    end

                end#end of if S_org2


            end #end of if seedN

        end #end of for ii

    end #end of loop N


    mx_idx = maximum(idx_coarse)

    fd_ssz = findall(x->x==0, idx_coarse)

    vec1 = collect(mx_idx + 1 : mx_idx + length(fd_ssz))

    idx_coarse[fd_ssz] = vec1

    push!(idx_mat, idx_coarse)

    ar_new = Any[]

    for i =1:length(ar)

        nd = ar[i]

        new_nd = idx_coarse[nd]

        push!(ar_new, new_nd)


    end #end of for i

    SZ = maximum(idx_coarse) + length(findall(x->x==0, idx_coarse))

    #append!(avg_cond, mean(Cvec_FB))
    append!(avg_cond, (sum(Cvec_FB)) / length(Cvec_FB))

    append!(sz_mat, SZ)

    append!(RR, (mx - SZ) / mx * 100)

    return avg_cond, ar_new

end #end of function

