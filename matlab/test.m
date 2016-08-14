%%%tuned to deal with rand_1000_10u,
%%newline issues if want to port back

clear

iters = 100;

for treeIndex = [0]
	for precision = [64,128,256,512,1024]
		digits(precision)

		folder = '../graphs/pathDisjoint_1000_exp20/';

		logName = strcat('_log_matlab_tree', num2str(treeIndex), '_', num2str(precision), 'precision_', num2str(iters), 'iters.txt')
		f_log = fopen(logName, 'w');

		fprintf(f_log, '=========CG using MATLAB variable precision=====\n');
		fprintf(f_log, '====NUMBER OF DIGITS = %d\n', digits);
		fprintf(f_log, '====DATA = %s\n', folder);
		fprintf(f_log, '====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);

		fprintf('=========CG using MATLAB variable precision=====\n');
		fprintf('====NUMBER OF DIGITS = %d\n', digits);
		fprintf('====DATA = %s\n', folder);
		fprintf('====error given as ||Ax - b||_2 / ||b||_2 and ||x - xbar||_A / ||xbar||_A\n', folder);

		f_matrix = strcat(folder, 'graph.mtx');
		f_tree = strcat(folder, 'tree', num2str(treeIndex), '.mtx');
		%f_tree = strcat(folder, 'tree.mtx');
		f_b = strcat(folder, 'b.vec');
		f_x = strcat(folder, 'x.vec');

		LG_default = getLaplacian(f_matrix);
		LG = vpa(LG_default);
		LT_default = getLaplacian(f_tree);
		LT = vpa(LT_default);

		n = size(LG, 1);

		A = full(LT_default < 0);
		T = graph(A);
		ord = bfsearch(T, 1);
		parent = zeros(n, 1);
		wParent = vpa(zeros(n, 1));

		seen = zeros(n, 1);
		seen(ord(1)) = 1;
		for i = 2:n
		  for j = 1:n
		      if seen(j) == 1 && A(ord(i), j) > 0
		          parent(i) = j;
		          wParent(i) = -vpa(LT_default(ord(i), j));
		      end
		  end
		  seen(ord(i)) = 1;
		end

		b = vpa(getVector(f_b));
		xbar = vpa(getVector(f_x));

		n = size(b, 1);
		onesN = vpa(ones(n, 1) / n);

		xbar = xbar - sum(xbar) * onesN; 
		sum(xbar)

		b = LG * xbar;

		b2 = norm(b);
		xbarA = sqrt(xbar' * LG * xbar);

		fprintf(f_log, '||b||_2 = %f\n', b2);
		fprintf(f_log, '||xbar||_A = %f\n', xbarA);
		fprintf('||b||_2 = %f\n', b2);
		fprintf('||xbar||_A = %f\n', xbarA);

		fclose(f_log);

		%[sum(xbar), sum(b)]

		%%%%CG copied from wiki
		x = vpa(zeros(n, 1));
		r = vpa(b);
		z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];
		z = z - sum(z) * onesN;
		p = z;

		for iter = 1:iters
		    alpha = (r' * z) / (p' * LG * p);
		    x = x + alpha * p;   
		   
		    rPrev = r;
		    zPrev = z;
		    
		    r = r - alpha * LG * p; %b - LG * x;
		    
		    f_log = fopen(logName, 'a');
		    fprintf(f_log,'i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, sqrt((x - xbar)'*LG*(x - xbar)) / xbarA);
		    fprintf('i=%3d, err2=%0.6g, errA = %0.6g\n', iter, norm(r) / b2, sqrt((x - xbar)'*LG*(x - xbar)) / xbarA);
		    fclose(f_log);
		    
		%    z = [LT(1:n-1,1:n-1) \ r(1:n-1); 0];

            z = r;
            %%Gauss backwards
            for i = n:-1:2
              z(parent(i)) = z(parent(i)) + z(ord(i));
            end
            %%Gauss down
            for i = 2:n
              z(ord(i)) = z(parent(i)) + z(ord(i)) / wParent(i);
            end
            z = z - sum(z) * onesN;   
		    
		    norm(LT * z - r)
		%[sum(z), sum(p)]
		    beta = (z' * r) / (zPrev' * rPrev);

		    p = z + beta * p;

		%[alpha, beta]
		end
    end

    
 
    norm(LT * z - r)

%[sum(z), sum(p)]
    beta = (z' * r) / (zPrev' * rPrev);

    p = z + beta * p;
end


