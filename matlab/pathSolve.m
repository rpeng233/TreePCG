%%%ARBITRARY PRECISION VERSION OF:
%z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
%z = z - sum(z) * onesN;
%assuming that 

tic

z = r;
for i = n:-1:2
    z(i - 1) = z(i - 1) + z(i);
end

for i = 2:n
    z(i) = (z(i - 1) - z(i)) / tree(i - 1, 3);
end
z = z - vpa(sum(vpa(z))) * onesN;

norm(LT * z - r)

toc