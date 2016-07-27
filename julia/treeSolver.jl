" An exact solver for a tree system. Takes in the adjacency matrix of the tree. "
function treeSolver{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti})
    
    n = tree.n;

    ord = dfsOrder(tree)
    permTree = tree[ord,ord]
    permLapTree = lap(tree[ord,ord])
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
        aux = copy(b[ord]) - mean(b);

        for u in n:-1:1
            for ind in 1:deg(permTree,u)
                v = nbri(permTree,u,ind)
                
                if v < u
                    fact = 1
                    geld[v] -= geld[u]
                    aux[v] += aux[u]
                end
            end
        end

        res = ones(Tval,n);
        for i in 2:n
            res[i] = (aux[i] + geld[i] * res[father[i]]) / geld[i]
        end
        res = res - mean(res);
        res = res[invperm(ord)];

        return res
    end

    return f

end