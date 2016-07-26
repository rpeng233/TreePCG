function toHighPrecision{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}; precision::Int64=1024)
    u,v,w = findnz(mat)
    w1024 = Array{BigFloat,1}(0)
    for val in w
        push!(w1024, BigFloat(val))
    end
    return sparse(u,v,w1024)
end

function toHighPrecision{Tv}(v::Array{Tv,1}; precision::Int64=1024)
    res = Array{BigFloat,1}(0)
    for val in v
        push!(res, BigFloat(val))
    end
    return res
end