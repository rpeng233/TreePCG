clear

digits(1024);

folder = '../graphs/pathDisjoint_30_exp20/';
f_log = fopen('log.txt', 'w');

fprintf(f_log, '=========CG using MATLAB variable precision=====\n');
fprintf(f_log, '====NUMBER OF DIGITS = %d\n', digits);
fprintf(f_log, '====DATA = %s\n', folder);
fprintf(f_log, '====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);

fprintf('=========CG using MATLAB variable precision=====\n');
fprintf('====NUMBER OF DIGITS = %d\n', digits);
fprintf('====DATA = %s\n', folder);
fprintf('====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);
fclose(f_log);

f_matrix = strcat(folder, 'graph.mtx');
f_tree = strcat(folder, 'tree.mtx');
f_b = strcat(folder, 'b.vec');
f_x = strcat(folder, 'x.vec');

LG = vpa(getLaplacian(f_matrix));
LT = vpa(getLaplacian(f_tree));

b = vpa(getVector(f_b));
xbar = vpa(getVector(f_x));
n = size(b, 1);
onesN = vpa(ones(n, 1) / n);

b2 = norm(b);
xbarA = sqrt(xbar' * LG * xbar);
fprintf('||b||_2 = %f\n', b2);
fprintf('||xbar||_A = %f\n', xbarA);

%[sum(xbar), sum(b)]

%%%%CG copied from wiki
x = vpa(zeros(n, 1));
r = vpa(b);
z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
z = z - sum(z) * onesN;
p = z;

for iter = 1:1000
    alpha = (r' * z) / (p' * LG * p);
    x = x + alpha * p;   
    
    rPrev = r;
    zPrev = z;
    
    r = r - alpha * LG * p; %b - LG * x;
    
    f_log = fopen('log.txt', 'a');
    fprintf(f_log,'i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, (x - xbar)'*LG*(x - xbar) / xbarA);
    fprintf('i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, (x - xbar)'*LG*(x - xbar) / xbarA);
    fclose(f_log);
    
    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
    z = z - sum(z) * onesN;
    
%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
%[alpha, beta]
end




