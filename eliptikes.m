%%
clear;
clc;
 
% bisn per dimension 
n=100;

% mesh size per dimension
h=1/n;

% indentity matrices
I=speye(n-1); % to speye δεν κρατα μηδενικα 

e=ones(n-1,1);
T=spdiags([4*e -10*e 4*e],[-1 0 1],n-1,n-1);
T2=spdiags([1*e 1*e], [-1 1],n-1,n-1);

%create the matrix A 
A2D=1/(6*h^2)*(kron(I,T)+kron(T,I)+kron(T2,T2));

% form rhs
f=inline('-2*pi^2*sin(pi*x).*sin(pi*y)');
%solve 
x=h:h:1-h;
y=h:h:1-h;

[xx,yy]=ndgrid(x,y);
F=f(xx(:),yy(:));

D=-gallery('poisson',n-1);
x=A2D\(F+1/12*D*F);

mesh(reshape(x,n-1,n-1));
figure
ex=sin(pi*xx(:)).*sin(pi*yy(:));
mesh(reshape(abs(x-ex),n-1,n-1));

