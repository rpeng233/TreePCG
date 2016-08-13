%%reads in a list of directed edges
%%  with weights stored as conductances
%%outputs an undirected graph Laplacian
function L = getLaplacian(fileName) 
    nonZeros = dlmread(fileName, ' ', 1, 0);
    n = nonZeros(1, 1);
    m = size(nonZeros, 1);
    A = sparse(nonZeros(2:m, 1), nonZeros(2:m, 2), nonZeros(2:m, 3), n, n, m - 1);
    A = A + A';
    D = sparse(1:n,1:n,sum(A)');
    L = D - A;
end
