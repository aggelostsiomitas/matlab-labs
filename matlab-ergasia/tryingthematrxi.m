clear;
clc;
n=25;
 
A=2*eye(n);
A(1:n-1,2:n)= A(1:n-1,2:n)-eye(n-1);
A(1:n-2,3:n)=A(1:n-2,3:n)-eye(n-2);
A(2:n,1:n-1)=A(2:n,1:n-1)-2*eye(n-1);
A(3:n,1:n-2)=A(3:n,1:n-2)-2*eye(n-2);
N=rand(n,1);
disp(A);
