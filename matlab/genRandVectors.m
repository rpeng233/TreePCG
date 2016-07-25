x = rand(n, 1);
x = x - mean(x) * ones(n, 1);
b = LG * x;

dlmwrite(f_b, double(b),'precision','%.10f');
dlmwrite(f_x, double(x),'precision', '%.10f');