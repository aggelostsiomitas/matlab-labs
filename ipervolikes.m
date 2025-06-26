%%
clear;
clc;

%number of intervals 
n=200;
nt=400;
c=0.13;

%mesh size
h=1/n;

%Laplacian
A=-gallery('tridiag',n-1); 
I=speye(n-1);

%time step 
Tmax=4;
dt=Tmax/nt;


%initial condition 
x=linspace(0,1,n+1)';

u=zeros(n+1,n+1);
u(:,1)=sin(3*pi*x);

% compute r
r=c*dt/h;

%first step for t=0
u(2:n,2)=0.5*(2*I+r^2*A)*u(2:n,1);

%the remaining steps 
for j=3:nt+1
    u(2:n,j)=(2*I+r^2*A)*u(2:n,j-1)-u(2:n,j-2);
end

%plot 
mesh(u);
%%


