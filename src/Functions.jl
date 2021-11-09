## Function to read inout and list the nodes belong to each hyperedge

function ReadInp(input)

    io = open(input, "r")
    ar  = Any[]

    while !eof(io)
        rr = zeros(Int, 0)
        ln = readline(io)
        sp = split(ln)

        for kk = 1:length(sp)
            r = parse(Int, sp[kk])
            append!(rr, r)
        end #kk
        push!(ar, rr)

    end

    ar = deleteat!(ar, 1)

    return ar

end #end of function


# list the hyperedges belong to each node
function HyperNodes(ar::Array{Any,1})

    H = INC3(ar)

    NH1 = Any[]

    rr1 = H.rowval

    cc1 = H.colptr

    for i = 1:size(H, 2)

        st = cc1[i]

        ed = cc1[i+1] - 1

        push!(NH1, rr1[st:ed])

    end

    return NH1

end


# create the incidence matrix of the input hypergraph
function INC3(ar::Array{Any,1})

    col = zeros(Int, 0)
    row = zeros(Int, 0)

    for iter = 1:length(ar)
        cc = (iter) * ones(Int, length(ar[iter]))
        rr = ar[iter]

        append!(col, cc)
        append!(row, rr)
    end

    row = row

    val = ones(Float64, length(row))

    mat = sparse(col, row, val)

    return mat
end

# find the maximum value in the array
function mx_func(ar::Array{Any,1})

    mx2 = Int(0)
    aa = Int(0)

    for i =1:length(ar)

    	mx2 = max(aa, maximum(ar[i]))
    	aa = mx2

    end

    return mx2

end

# convert the hypergraph to asimple graph applying Star expansion
function Star(ar::Array{Any,1})

    mx = mx_func(ar)

    sz = length(ar)

    col = zeros(Int32, 0)
    val = zeros(Float32, 0)
    row = zeros(Int32, 0)

    for iter =1:length(ar)
        LN = length(ar[iter])
        cc = (iter+mx) * ones(Int, LN)
        vv = (1/LN) * ones(Int, LN)

        rr = ar[iter]
        append!(col, cc)

        append!(row, rr)

        append!(val, vv)
    end

    mat = sparse(row, col, val,mx+sz, mx+sz)

    A = mat + mat'

    return A

end

# Filter high frequency components of the input vectors
function Filter(sm_vec::Array{Float64,2}, k::Int64, AD::SparseMatrixCSC{Float32,Int32}, mx::Int64)

    N = size(sm_vec, 2);

    sz = size(AD, 1)

    V = zeros(mx, N);

    AD = AD .* 1.0

    AD[diagind(AD, 0)] = AD[diagind(AD, 0)] .+ 0.1

    dg = sum(AD, dims = 1) .^ (-.5)

    I2 = 1:sz

    D = sparse(I2, I2, sparsevec(dg))

    for iter in 1:N

        sm = sm_vec[:, iter]

        for loop in 1:k

            sm = D * sm

            sm = AD * sm

            sm = D * sm

        end

        sm = sm[1:mx]

        on = ones(Int, mx)

        #Make all the eigenvectors orthogonal to all 1 vectors
        sm_ot = sm - ((dot(sm, on) / dot(on, on)) * on)

        # Normalization
        sm_norm = sm_ot ./ norm(sm_ot);

        V[: , iter] = sm_norm;

    end

    AD[diagind(AD, 0)] .= 0
    AD = dropzeros!(AD)

    return V

end #end of function

# compute the hyperedge scores using non-linear quadratic form
# max function (without mediators)
function h_score3(ar::Array{Any,1}, SV::Array{Float64,2})
    score = zeros(eltype(SV), length(ar))
    @inbounds Threads.@threads for i in eachindex(ar)
        nodes = ar[i]
        for j in axes(SV, 2)
            mx, mn = -Inf, +Inf
            for node in nodes
                x = SV[node, j]
                mx = ifelse(x > mx, x, mx)
                mn = ifelse(x < mn, x, mn)
            end
            score[i] += (mx - mn)^2
        end
    end
    return score
end


# k-mean clustering inside each hyperedge using
#the spectral simple graph (star graph) embeddings
function he_cluster14(idx::Array{Int64,1}, ar::Array{Any,1}, Hscore::Array{Float64,1}, SV::Array{Float64,2}, mx::Int64)

    ratio = 2

    val = 1

    MX = maximum(Hscore);

    M = length(ar)

    EL = collect(UnitRange(1,M))

    idx_new = zeros(Int, mx)

    new_NL = Any[]

    new_NN = Any[]

    rem_vec = []

    flag = falses(1, mx)

    FE = falses(M)

    ## node clustering

    ss = zeros(Int, 0)

    for i = 1:length(ar)

        append!(ss, length(ar[i]))

    end

    E = sortperm(ss, rev=true)

    for iter2 = 1:length(E)

        nd = ar[E[iter2]]

        fd_flag = flag[nd]

        node = nd[.!fd_flag]


        if length(node) > 0

            sm_mat = view(SV, node, :)

            R = kmeans(sm_mat', ceil(Int, length(node)/ratio));

            idx_clus = R.assignments

            idx_new[node] = idx_clus .+ val;

            flag[node] .= 1

            val = val + maximum(idx_clus);

        end # end if length(node) > 1


    end # end of iter2

    return idx_new

end#function

# finding non-unique clusters
function nonunique2!(x::AbstractArray{T}) where T
    x = sort(x)
    duplicatedvector = T[]
    for i=2:length(x)
        if (isequal(x[i],x[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], x[i])))
            push!(duplicatedvector,x[i])
        end
    end
    duplicatedvector
end
  
    
# Write the output matrix in hMetis format
function Whgr(input, ar)
    mx = mx_func(ar)
    open(input,"w")do io
        println(io, length(ar)," ", mx)
        for i =1:length(ar)
            nds = ar[i]
            for j =1:length(nds)
                print(io, nds[j], " ")
            end
            println(io)
        end
    end
end
