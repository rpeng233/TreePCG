%%reads in a list of directed edges
%%  with weights stored as conductances
%%outputs an undirected graph Laplacian
function edges = getEdges(fileName) 
    edges = dlmread(fileName, ' ', 2, 0);
end