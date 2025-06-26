clc; 
clear;

t1 = tic; %start time countdown

% Grid size
m = 8;

% Domain
L = [0, 0.4];

% Step size
h = (L(2) - L(1)) / m; % 0.05
hy=h;
hx=h;

% Number of positions
Nx = m - 1;
Ny = m - 1;
n = Nx * Ny;

% Grid points
x = linspace(L(1) + hx, L(2) - hx, Nx);
y = linspace(L(1) + hy, L(2) - hy, Ny);
[xx,yy]=ndgrid(x,y);
xx = reshape(xx, n, 1);
yy = reshape(yy, n, 1);

% matrixes with ones
Ix = speye(Nx);
Iy = speye(Ny);

% Tridiagonal matrixes
Tx = -gallery('tridiag', Nx) / hx^2;
Ty = -gallery('tridiag', Ny) / hy^2;

% matrix A
A = kron(Iy, Tx) + kron(Ty, Ix) + (12.5 * pi^2 *h^2* speye(n)); 

% f function
f=@(x,y) -25*pi^2*sin((5*pi/2).*x).*sin((5*pi/2).*y);
F = f(xx, yy);

% Solve for numerical solution
x_sol = A \ F;

% Plot numerical solution
figure;
mesh(reshape(x_sol, Nx, Ny));
title('Numerical Solution');

% Analytical solution
U_ANAL = sin(5 * pi .* xx / 2) .* sin(5 * pi .* yy / 2);

% Absolute error plot
figure;
mesh(reshape(abs(x_sol - U_ANAL), Nx, Ny));
title('Absolute Error');


% 9-Point Method 
e = ones(Nx, 1);

T1 = spdiags([4 * e, (-10 - 3 * h^2 * 12.5 * pi^2) * e, 4 * e], [-1, 0, 1], Nx, Ny);
T2=spdiags([e,e],[-1,1],Nx,Ny);
A2D=1/(6*h^2)*(kron(Ix,T1)+kron(T1,Ix)+kron(T2,T2));

D = -gallery('poisson', Nx);
F = f(xx, yy);

%solve
x2 = A2D \ (F + 1/12 * D * F);

% Plot 9-Point Method numerical solution
figure;
mesh(reshape(x2, Nx, Ny));
title('Numerical Solution - 9-Point Method');

% Absolute error for 9-Point Method
figure;
mesh(reshape(abs(x2 - U_ANAL), Nx, Ny));
title('Absolute Error - 9-Point Method');

%end time
t2 = toc(t1);
disp(['Συνολικός χρόνος που χρειάστηκε για τους υπολογισμούς: ', num2str(t2)]);

