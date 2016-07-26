""" 
    Generates a chain with edges between i and i+2, if i equals 1 modulo k. 
    Returns both the graph and the chain as its corresponding tree.
"""
function chainCycle(n::Int64; k::Int64=3)
    u = Array{Int64,1}(0)
    v = Array{Int64,1}(0)
    w = Array{Float64,1}(0)

    for i in 1:(n-1)
        push!(u, i)
        push!(v, i+1)
        push!(w, 1.0)
        
        push!(u, i+1)
        push!(v, i)
        push!(w, 1.0)
    end

    valpool = rand() * exp(20 * rand(n));

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

    a = sparse(u,v,w); la = lap(a); n = a.n;

    utree = collect(1:(n-1))
    append!(utree, collect(2:n))

    vtree = collect(2:n)
    append!(vtree, collect(1:(n-1)))

    wtree = ones(2 * n - 2)
    tree = sparse(utree,vtree,wtree)

    return a, tree;

end