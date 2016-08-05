" An exact solver for a tree system. Takes in the adjacency matrix of the tree. "
function treeSolver{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}; d::Array{Tv,1} = zeros(Tv,tree.n))
    
    n = tree.n;

    ord = dfsOrder(tree)
    permTree = tree[ord,ord]
    permLapTree = lap(permTree) + spdiagm(d[ord])
    father = ones(Int64,n);

    for u in 2:n
        for i in 1:deg(permTree,u)
            v = nbri(permTree,u,i)
            if v < u
                father[u] = v
                break
            end
        end
    end

    function f{Tval}(b::Array{Tval,1})

        geld = diag(permLapTree)
        aux = copy(b[ord])
        if norm(d) == 0
            aux = aux - mean(b);
        end

        for u in n:-1:1
            for ind in 1:deg(permTree,u)
                v = nbri(permTree,u,ind)
                w = weighti(permTree,u,ind)
                
                if v < u
                    fact = w / geld[u]

                    geld[v] -= w * fact
                    aux[v] += aux[u] * fact
                end
            end
        end

        res = ones(Tval,n);
        if norm(d) != 0
            res[1] = aux[1] / geld[1]
        end

        for i in 2:n
            res[i] = (aux[i] + permTree[father[i],i] * res[father[i]]) / geld[i]
        end

        if norm(d) == 0
            res = res - mean(res)
        end
        res = res[invperm(ord)];

        return res
    end

    return f

end