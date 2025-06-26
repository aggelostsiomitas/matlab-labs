clear;
clc;
n =1000;

A = 2*eye(n);
A(1:n-1,2:n)=A(1:n-1,2:n)-eye(n-1); 
A(1:n-2,3:n)=A(1:n-2,3:n)-eye(n-2);
A(2:n,1:n-1)=A(2:n,1:n-1)-eye(n-1);
A(3:n,1:n-2)=A(3:n,1:n-2)-eye(n-2);

N = rand(n, 1);
figure;
disp(A);
spy(A);

[c_result, elapsed_time] = multiply_A(A, N,n);
figure;
disp('Result of multiply_A:');
disp(c_result);
spy(c_result);
fprintf('Time elapsed for multiply_A: %.6f seconds\n', elapsed_time);

[x_result, elapsed_time2] = multiply_AT(A,N,n);
figure;
disp('Result of multiply_AT:');
disp(x_result);
spy(x_result);
fprintf('Time elapsed for multiply_AT: %.6f seconds\n', elapsed_time2);
