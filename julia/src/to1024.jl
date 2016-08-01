function to1024{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
    u,v,w = findnz(mat)
    w1024 = Array{BigFloat,1}(0)
    for val in w
        push!(w1024, BigFloat(val))
    end
    return sparse(u,v,w1024)
end

function to1024{Tv}(v::Array{Tv,1})
    res = Array{BigFloat,1}(0)
    for val in v
        push!(res, BigFloat(val))
    end
    return res
end