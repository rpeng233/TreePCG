function randGraph(n::Int64; p::Float64=0.5, weightGen::Function=rand)
    u = Array{Int64,1}(0)
    v = Array{Int64,1}(0)
    w = Array{Float64,1}(0)

    nr = ceil(Int64, n * p / 2)
    for i in 1:n
        samps = rand(1:n, nr)
        iv = Array{Int64,1}(0)
        for j in samps
            if j != i
                push!(iv, j)
            end
        end

        iu = [i for j in 1:length(iv)]
        iw = ones(length(iv))

        append!(u, iu)
        append!(v, iv)
        append!(w, iw)
    end
    mat = sparse(u,v,w)

    if mat.m != mat.n
        println("Failed - graph is disconnected")
        return -1
    end
    mat = mat + mat'

    # get proper weights
    u,v,w = findnz(tril(mat))
    w = Float64[weightGen() for i in 1:length(w)]
    mat = sparse(u,v,w,n,n);
    mat = mat + mat'

    if maximum(components(mat)) > 1
        println("Failed - graph is disconnected")
        return -1
    end

    return mat
end


""" 
    Generates a chain with edges between i and i+2, if i equals 1 modulo k. 
    Weights on the chain are ~Unif(0,1). Weights on the cycle edges are ~Unif(0,1) * e^(30 * Unif(0,1))
    Returns both the graph and the chain as its corresponding tree.
"""
function chainCycle(n::Int64; k::Int64=3)
    u = Array{Int64,1}(0)
    v = Array{Int64,1}(0)
    w = Array{Float64,1}(0)

    for i in 1:(n-1)
        val = rand()
        push!(u, i)
        push!(v, i+1)
        push!(w, val)
        
        push!(u, i+1)
        push!(v, i)
        push!(w, val)
    end

    utree = copy(u)
    vtree = copy(v)
    wtree = copy(w)
    tree = sparse(utree,vtree,wtree)

    valpool = rand() * exp(30 * rand(n));

    for i in 1:(n-2)
        if i % k == 1
            val = valpool[i]
            
            push!(u, i)
            push!(v, i+2)
            push!(w, val)
            
            push!(u, i+2)
            push!(v, i)
            push!(w, val)
        end
    end

    a = sparse(u,v,w);

    return a, tree;

end
