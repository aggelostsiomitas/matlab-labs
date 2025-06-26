0%%
clear;
clc;
 
% bisn per dimension 
nx=5;
ny=5;


% mesh size per dimension
hx=1/nx;
hy=1/ny;

% indentity matrices
Ix=speye(nx-1); % to speye δεν κρατα μηδενικα 
Iy=speye(ny-1);

Tx=-gallery('tridiag',nx-1)/hx^2;
Ty=-gallery('tridiag',ny-1)/hy^2;

%create the matrix A 
A2D=kron(Iy,Tx)+kron(Ty,Ix);

% form rhs
F=ones((nx-1)*(ny-1),1);

%solve 
x=A2D\F;
    
mesh(reshape(x,nx-1,ny-1));

%%
clear;
clc;
 
% bisn per dimension 
nx=20;
ny=20;
nz=20;

% mesh size per dimension
hx=1/nx;
hy=1/ny;
hz=1/nz;
% indentity matrices
Ix=speye(nx-1); % to speye δεν κρατα μηδενικα 
Iy=speye(ny-1);
Iz=speye(nz-1);

Tx=-gallery('tridiag',nx-1)/hx^2;
Ty=-gallery('tridiag',ny-1)/hy^2;
Tz=-gallery('tridiag',nz-1)/hy^2;
%create the matrix A 
A3D=kron(Iz,kron(Iy,Tx))+ kron(Iz,kron(Ty,Ix))+kron(Tz,kron(Iy,Ix));


% form rhs
F=ones((nx-1)*(ny-1)*(nz-1),1);

%solve 
x=A3D\F;
    
xr=reshape(x,nx-1,ny-1,nz-1);

for i=1:nz-1
    mesh(xr(:,:,i));
    M=getframe;
    pause(0.3);
end
% an eixa d^2u/dx^2+d^u/dy^2+du/dx=1 
%to du/dx exei morio  [-1 0 1]

%%
clear;
clc;
 
% bisn per dimension 
nx=50;
ny=50;


% mesh size per dimension
hx=1/nx;
hy=1/ny;

% indentity matrices
Ix=speye(nx); % to speye δεν κρατα μηδενικα 
Iy=speye(ny-1);

Tx=-gallery('tridiag',nx)/hx^2;
Tx(nx,nx-1)=2/hx^2;
Ty=-gallery('tridiag',ny-1)/hy^2;

%create the matrix A 
A2D=kron(Iy,Tx)+kron(Ty,Ix);

% form rhs
F=ones((nx)*(ny-1),1);

%solve 
x=A2D\F;
    
mesh(reshape(x,nx,ny-1));
%%
clear;
clc;
 
% bisn per dimension 
nx=50;
ny=50;


% mesh size per dimension
hx=1/nx;
hy=1/ny;

% indentity matrices
Ix=speye(nx); % to speye δεν κρατα μηδενικα 
Iy=speye(ny);

Tx=-gallery('tridiag',nx)/hx^2;
Tx(nx,nx-1)=2/hx^2;
Ty=-gallery('tridiag',ny)/hy^2;
Ty(ny,ny-1)=2/hy^2;
%create the matrix A 
A2D=kron(Iy,Tx)+kron(Ty,Ix);

% form rhs
F=ones((nx)*(ny),1);

%solve 
x=A2D\F;
    
mesh(reshape(x,nx,ny));