clear

f_log = fopen('log.txt', 'w');

folder = '../graphs/se_2/';

fprintf(f_log, '=========CG using MATLAB default precision=====\n');
fprintf(f_log, '====DATA = %s\n', folder);
fprintf(f_log, '====error given as ||Ax - b||_2^2 and ||x - xbar||_A^2\n', folder);

fprintf('=========CG using MATLAB default precision=====\n');
fprintf('====DATA = %s\n', folder);
fprintf('====error given as ||Ax - b||_2^2 and ||x - xbar||_A^2\n', folder);

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

for iter = 1:10000
    alpha = (r' * z) / (p' * LG * p);
    x = x + alpha * p;
   
    
    rPrev = r;
    zPrev = z;
    
    r = r - alpha * LG * p; %b - LG * x;

    
    fprintf(f_log, 'i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r)^2, (x - xbar)'*LG*(x - xbar));
    fprintf('i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r)^2, (x - xbar)'*LG*(x - xbar));
    
    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
    z = z - sum(z) * onesN;
%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
%[alpha, beta]
end

fclose(f_log);

