function [x, elapsed_time2] = multiply_AT(A, N, n)
  D= 3; 
  x = zeros(n, 1);
  t2 = tic;
  for i = 1:n
      for j = max(i-D+1, 1):min(i+D-1, n)
          x(i) = x(i) + A(i, j)' * N(j); 
      end
  end
  
  elapsed_time2 = toc(t2); 
  disp(elapsed_time2);
end
