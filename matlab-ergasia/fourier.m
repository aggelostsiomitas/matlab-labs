clear;
clc;
n=1:1:300;
pi=3.14;
a0=0;
A=@(n)(1/(pi*n))*sin(0.8*pi*n)+(1/(0.8*pi^2*n^2))*(cos(0.8*pi*n)-2+cos(1.2*pi*n));
B=@(n)(1/(pi*n))*(-cos(0.8*pi*n)+1) +(1/(0.8*pi^2*n^2))*(sin(0.8*pi*n)-sin(1.2*pi*n));
A_vals = zeros(1, length(n));
B_vals = zeros(1, length(n));
C_vals = zeros(1, length(n));
A_vals(1) = 0; 
C_vals(1) = A_vals(1); 
for i = 2:length(n)
    A_vals(i) = A(i);
    B_vals(i) = B(i);
    C_vals(i) = sqrt(A_vals(i).^2 + B_vals(i).^2);
end


hold on;
plot(n,A_vals,'b');
plot(n,B_vals,'r');
plot(n,C_vals,'g');
hold off;