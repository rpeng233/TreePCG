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

        if verbose && itcnt % 10 == 0
            println("Working on iteration ", itcnt, " with a-norm error ", Float64(sqrt(abs((lhs - x)' * mat * (lhs - x))[1])))
        end
        
        q = mat*p

        al = rho/dot(p, q)

        x = x + al * p
        r -= al * q
        # BLAS.axpy!(al,p,x)  # x = x + al * p
        # BLAS.axpy!(-al,q,r)  # r -= al*q

        if true
            errMN = Float64(sqrt(abs((lhs - x)' * mat * (lhs - x))[1]))
            err2_a = Float64(norm(mat * x - b) / norm(b))
            err2_b = Float64(norm(lhs - x))
            push!(dbg, "iter = $(itcnt)      errA=$(errMN)      err2_a=$(err2_a)      err2_b=$(err2_b)")
        end

        # here is the top of the code in numerical templates
        z = pre(r)

        oldrho = rho
        rho = dot(z, r)

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho

        p = z + beta*p
        # BLAS.scal!(n,beta,p,1) # p *= beta
        # BLAS.axpy!(1.0,z,p) # p += z
      end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", Float64(sqrt(abs((lhs - x)' * mat * (lhs - x))[1])), ".")
    end

    return x,dbg
end
