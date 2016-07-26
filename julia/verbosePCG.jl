function pcgV{Tv}(mat, b::Array{Tv,1}, pre::Function; tol::Tv=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    pcgVBLAS(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

function pcgVMatNorm{Tv}(mat, b::Array{Tv,1}, pre::Function, lhs::Array{Tv,1}; tol::Tv=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    return pcgVBLASMatNorm(mat, b, pre, lhs, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function pcgVBLAS{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Tval=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1

        if verbose && itcnt % 100 == 0
            println("Working on iteration ", itcnt)
        end
        
        q = mat*p

        al = rho/dot(p, q)

        x = x + al * p
        # BLAS.axpy!(al,p,x)  # x = x + al * p

        r = r - al * q
        # BLAS.axpy!(-al,q,r)  # r -= al*q

        if (verbose && (itcnt < 100 || itcnt % 100 == 0)) ||
           (norm(r) < tol)
            err = Float64(norm(r) / norm(b))
            push!(dbg, "$(err) error after $(itcnt) iterations")
        end

        if norm(r) < tol #Converged?
            break
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

        if (time() - t1) > maxtime
            if verbose
                println("PCG stopped at maxtime.")
            end
            break
        end

       
      end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    return x
end


function pcgVBLASMatNorm{Tval}(mat, b::Array{Tval,1}, pre, lhs::Array{Tval,1}; tol::Tval=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1

        if verbose && itcnt % 100 == 0
            println("Working on iteration ", itcnt, " with error ", sqrt((lhs - x)' * mat * (lhs - x))[1] )
        end
        
        q = mat*p

        al = rho/dot(p, q)

        x = x + al * p
        r -= al * q
        # BLAS.axpy!(al,p,x)  # x = x + al * p
        # BLAS.axpy!(-al,q,r)  # r -= al*q

        if (verbose && (itcnt < 100 || itcnt % 100 == 0)) ||
           (sqrt((lhs - x)' * mat * (lhs - x))[1] < tol)
            err = Float64(sqrt((lhs - x)' * mat * (lhs - x))[1])
            push!(dbg, "$(err) error(matrix norm) after $(itcnt) iterations")
        end

        # This will take some time. Can we do something faster?
        if sqrt((lhs - x)' * mat * (lhs - x))[1] < tol 
            break
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

        if (time() - t1) > maxtime
            if verbose
                println("PCG stopped at maxtime.")
            end
            break
        end

      end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", sqrt((lhs - x)' * mat * (lhs - x))[1], ".")
    end

    return x
end
