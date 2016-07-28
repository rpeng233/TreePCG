""" 
    Generates a chain with edges between i and i+2, if i equals 1 modulo k. 
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

    valpool = rand() * exp(40 * rand(n));

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