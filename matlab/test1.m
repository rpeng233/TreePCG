clear

digits(10);

folder = '../graphs/se/1/';

f_matrix = strcat(folder, 'graph.mtx');
f_tree = strcat(folder, 'tree.mtx');
f_vec = strcat(folder, 'b.vec');
f_x = strcat(folder, 'x.vec');

graph = getEdges(f_matrix);
LG = getLaplacian(f_matrix);

tree = getEdges(f_tree);
LT = getLaplacian(f_tree);

b = vpa(getVector(f_vec));
xbar = vpa(getVector(f_x));
n = size(b, 1);
onesN = vpa(ones(n, 1) / n);

[sum(xbar), sum(b)]


%%%%CG copied from wiki
x = vpa(zeros(n, 1));
r = vpa(b);
pathSolve;
p = z;
% % 
% % for iter = 1:40
% %     alpha = (r' * z) / (p' * LG * p);
% %     x = x + alpha * p;   
% %     
% %     rPrev = r;
% %     zPrev = z;
% %     
% %     r = r - alpha * LG * p; %b - LG * x;
% % fprintf('iteration %3d, residue = %0.6g\n', iter, norm(r));
% %     
% %     pathSolve;
% % [sum(z), sum(p)]
% %     beta = (z' * r) / (zPrev' * rPrev);
% % 
% %     p = z + beta * p;
% % [alpha, beta]
% % end
% % 



