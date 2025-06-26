clear;
clc;

%dips of matrix 
n=10;

%create the linear system
A=gallery('tridiag',n);

%create the right side 
x_exact=(1:n)'/n;
b=A*x_exact;

%maximum number allowed
nmax=1e5;

%tolerance
tol=1e-5;

%extract the diagonal part of A 
M=diag(diag(A))^-1;

% to ena diag kanei mia stili kai to deutero to kanei se pinaka san
% kanonikh diagonio 

%initial guess
x0=zeros(n,1);

%absolute kritirio of termination 
normb=norm(b);


%start iterations
for i=1:nmax
   x=x0+M*(b-A*x0);
   if norm(b-A*x)<tol*normb %best tactic to drop the tol
       break;
   end
   x0=x;
end
if (i==nmax)
    disp(['No convergence after :',num2str(namx)]);
else
    disp(['iterations for  convergence :',num2str(i)]);
end

