function [c, elapsed_time] = multiply_A(A,N,n)
    D = 3;
    c = zeros(n, 1);
    t1 = tic; 
    for i = 1:n
        for j = max(i-D+1, 1):min(i+D-1, n)
            c(i) = c(i) + A(i, j) * N(j); 
        end
    end
    elapsed_time = toc(t1); 
    disp(elapsed_time);
end


