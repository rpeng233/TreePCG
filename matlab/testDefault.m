clear

fprintf('CG using MATLAB default precision\n');


folder = '../graphs/se/1/';
f_matrix = strcat(folder, 'graph.mtx');
f_tree = strcat(folder, 'tree.mtx');
f_b = strcat(folder, 'b.vec');
f_x = strcat(folder, 'x.vec');

LG = getLaplacian(f_matrix);
LT = getLaplacian(f_tree);
b = getVector(f_b);
xbar = getVector(f_x);

%[sum(xbar), sum(b), norm(LG * xbar - b)]

n = size(b, 1);
onesN = ones(n, 1) / n;

%%%%CG copied from wiki
x = zeros(n, 1);
r = b;
z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
z = z - sum(z) * onesN;    
p = z;

%z' * z

%p' * LG * p

for iter = 1:100
    alpha = (r' * z) / (p' * LG * p);
    x = x + alpha * p;
   
    
    rPrev = r;
    zPrev = z;
    
    r = r - alpha * LG * p; %b - LG * x;
fprintf('iteration %3d, residue = %0.6g\n', iter, norm(r));
    
    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
    z = z - sum(z) * onesN;
%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
%[alpha, beta]
end




