clear

digits(1024);

folder = '../graphs/se/1/';
f_matrix = strcat(folder, 'graph.mtx');
f_tree = strcat(folder, 'tree.mtx');
f_b = strcat(folder, 'b.vec');
f_x = strcat(folder, 'x.vec');

LG = getLaplacian(f_matrix);
LT = getLaplacian(f_tree);

b = vpa(getVector(f_b));
xbar = vpa(getVector(f_x));
n = size(b, 1);
onesN = vpa(ones(n, 1) / n);

%[sum(xbar), sum(b)]

digits

%%%%CG copied from wiki
x = vpa(zeros(n, 1));
r = vpa(b);
z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
z = z - sum(z) * onesN;
p = z;

for iter = 1:100
    alpha = (r' * z) / (p' * LG * p);
    x = x + alpha * p;   
    
    rPrev = r;
    zPrev = z;
    
    r = r - alpha * LG * p; %b - LG * x;
fprintf('iteration %3d, residue = %0.6g\n', iter, double(norm(r)));
    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
    z = z - sum(z) * onesN;
    
%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
%[alpha, beta]
end




