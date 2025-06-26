%%
clc;
clear;

%parameters
c=0.5;
nx=10;
nt=10;
th=0.5;
Tmax=5*0.42;
dt=Tmax/nt;
L=1;
h=L/nx;

%matrixes
I=speye(nx-1);
D=-gallery('tridiag',nx-1)/h^2;

%x step
x=linspace(0,L,nx+1)';

%the matrix u that i am searching for
U=zeros(nx+1,nt+1);
U(:,1)=sin(pi*x);
for i=2:nt+1
    U(2:nx,i)=(I-th*dt*c*D)\(I+(1-th)*dt*c*D)*U(2:nx,i-1);
    %this was created from nicholson equation for th =0.5
end 
mesh(U)


