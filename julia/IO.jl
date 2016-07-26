""" Writes a matrix to a file in MatrixMarket format """
function writeToFile{Tv,Ti}(filename::ASCIIString, a::SparseMatrixCSC{Tv,Ti})
    f = open(filename, "w")

    println(f, "%%MatrixMarket matrix coordinate real symmetric")
    
    println(f, a.n, " ", a.m, " ", ceil(Int64,length(a.nzval) / 2))

    mat = tril(a);

    pos = 1
    for i in 1:length(mat.nzval)
        while mat.colptr[pos + 1] <= i
            pos = pos + 1
        end
        
        println(f, mat.rowval[i], " ", pos, " ", mat.nzval[i])
    end
    
    close(f)
end


""" Writes a vector to a file in MatrixMarket format """
function writeToFile{Tv}(filename::ASCIIString, a::Array{Tv,1})
    f = open(filename, "w")
    
    if typeof(a[1]) != ASCIIString
        println(f, length(a))

        for i in 1:length(a)
            println(f, a[i])
        end
    else
        for i in 1:length(a)
            println(f, a[i])
        end
    end
    
    close(f)
end


""" Read a matrix from a file in MatrixMarket format """
function readFromFile(filename::ASCIIString)

    f = open(filename)

    lines = readlines(f)
   
    ismat = split(lines[1], ' ')
    if ismat[3] == "coordinate"
        ismat = true
    else
        ismat = false
    end

    start = 1
    for i in 1:length(lines)
        if lines[i][1] == '%'
            start = i + 1
        else
            break
        end
    end

    if ismat
        # We are reading a Sparse Matrix
        m = parse(Int64, split(lines[start], ' ')[1])
        n = parse(Int64, split(lines[start], ' ')[2])
        nonzrs = parse(Int64, split(lines[start], ' ')[3])

        u = Array{Int64,1}(0)
        v = Array{Int64,1}(0)
        w = Array{Float64,1}(0)

        for i in 1:nonzrs
            p = parse(Int64, split(lines[i+start], ' ')[1])
            q = parse(Int64, split(lines[i+start], ' ')[2])
            r = parse(Float64, split(lines[i+start], ' ')[3])

            push!(u, p)
            push!(v, q)
            push!(w, r)

            push!(u, q)
            push!(v, p)
            push!(w, r)
        end

        return sparse(u,v,w,m,n)
    else
        # We are reading a vector
        n = parse(Int64, split(lines[start], ' ')[1])
        u = Array{Float64,1}(0)
        for i in 1:n
            push!(u, parse(Float64, lines[i+start]))
        end

        return u
    end

    close(f)

end
