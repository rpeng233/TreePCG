clear

f_log = fopen('log.txt', 'w');

folder = '../graphs/pathDisjoint_1000_exp20/';

fprintf(f_log, '=========CG using MATLAB variable precision=====\n');
fprintf(f_log, '====NUMBER OF DIGITS = %d\n', digits);
fprintf(f_log, '====DATA = %s\n', folder);
fprintf(f_log, '====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);

fprintf('=========CG using MATLAB variable precision=====\n');
fprintf('====NUMBER OF DIGITS = %d\n', digits);
fprintf('====DATA = %s\n', folder);
fprintf('====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);

f_matrix = strcat(folder, 'graph.mtx');
f_tree = strcat(folder, 'tree.mtx');
f_b = strcat(folder, 'b.vec');
f_x = strcat(folder, 'x.vec');

LG = getLaplacian(f_matrix);
LT = getLaplacian(f_tree);
b = getVector(f_b);
xbar = getVector(f_x);

b2 = norm(b);
xbarA = sqrt(xbar' * LG * xbar);

fprintf(f_log, '||b||_2 = %f\n', b2);
fprintf(f_log, '||xbar||_A = %f\n', xbarA);
fprintf('||b||_2 = %f\n', b2);
fprintf('||xbar||_A = %f\n', xbarA);


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

for iter = 1:1000
    alpha = (r' * z) / (p' * LG * p);
    x = x + alpha * p;
   
    
    rPrev = r;
    zPrev = z;
    
    r = r - alpha * LG * p; %b - LG * x;

    fprintf(f_log,'i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, sqrt((x - xbar)'*LG*(x - xbar)) / xbarA);
    fprintf('i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, sqrt((x - xbar)'*LG*(x - xbar)) / xbarA);
    
    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
    z = z - sum(z) * onesN;
%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
%[alpha, beta]
end

fclose(f_log);

