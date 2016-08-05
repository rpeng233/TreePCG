function pcgV{Tv}(mat, b::Array{Tv,1}, pre::Function, lhs::Array{Tv,1}; maxits=1000, verbose=false)
    return pcgWorker(mat, b, pre, lhs, maxits=maxits, verbose=verbose)
end

# runs pcg for a given number of iterations. outputs debug to a dbg[] array (bad, but works for now). debug consists of both 2norm and anorm
function pcgWorker{Tv}(mat, b::Array{Tv,1}, pre, lhs::Array{Tv,1}; maxits=1000, verbose=false)

    dbg = []

    n = size(mat,2)
    
    x = zeros(Tv,n)
    
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
        
        q = mat*p

        al = rho/dot(p, q)

        x = x + al * p
        r -= al * q

        errMN = Float64(sqrt(((lhs - x)' * mat * (lhs - x))[1] / (lhs' * mat * lhs)[1]))
        err2 = Float64(norm(mat * x - b) / norm(b))

        push!(dbg, "iter = $(itcnt)      errA=$(errMN)      err2=$(err2)")

        if verbose && (itcnt % 10 == 1 || itcnt == maxits)
            println("Finished iteration ", itcnt, " with errors ", errMN, " ", err2)
        end

        # here is the top of the code in numerical templates
        z = pre(r)

        oldrho = rho
        rho = dot(z, r)

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho

        p = z + beta*p
      end

    return x,dbg
end
