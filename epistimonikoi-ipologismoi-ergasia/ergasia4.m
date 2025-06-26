clear;
clc;

t1=tic;
% Grid size
m = 4;

% Number of positions
n = (m - 1)^2;

% Domain
L = [0, 0.4];

% Step size
h = (L(2) - L(1)) / (m-1);

%grid parameters
Nx = m - 1; Ny = m - 1; 


x = linspace(L(1) + h, L(2) - h, Nx);
y = linspace(L(1) + h, L(2) - h, Ny);
[xx, yy] = ndgrid(x, y);

xx = reshape(xx, n, 1);
yy = reshape(yy, n, 1);

% 5-Point Method
A = zeros(n, n);
for i = 1:n
    A(i,i) = 4 + 12.5 * pi^2 * h^2;  
end

% Off-diagonal elements
for i = 1:n-1
    A(i, i+1) = -1;
    A(i+1, i) = -1;
end

% Handling boundary conditions
for i = 1:Nx-1
    A(i*Nx, i*Nx + 1) = 0;
    A(i*Nx + 1, i*Nx) = 0;
end

for i = 1:n-Nx
    A(i, i+Nx) = -1;
    A(i+Nx, i) = -1;
end

f = @(x,y) -25 * pi^2 * sin(5 * pi/2) .* x .* sin(5 * pi/2) .* y;
b = h^2 * f(xx, yy);  


x = A \ b;

% Plot numerical solution
figure;
mesh(reshape(xx, Nx, Ny), reshape(yy, Nx, Ny), reshape(x, Nx, Ny));
title('Numerical Solution');

% Analytical solution
U_ANAL = sin(5 * pi/2) * xx .* sin(5 * pi/2) .* yy;

% Error plot
figure;
mesh(reshape(xx, Nx, Ny), reshape(yy, Nx, Ny), reshape(abs(x - U_ANAL), Nx, Ny));
title('Absolute Error');



% Μέθοδος 9 σημείων
I = speye(Nx); 
e = ones(Nx,1);

T1 = spdiags([2/3*e (-10/3 - 6*h^2*12.5*pi^2)*e 2/3*e], [-1 0 1], Nx, Nx);
T2 = spdiags([1*e 1*e], [-1 1], Nx, Nx);

A2D = 1/(6*h^2)*(kron(I,T1) + kron(T1,I) + kron(T2,T2));
F = f(xx, yy);

D=-gallery('poisson',Nx);
x2 = A2D \ (F + 1/12 * D * F);

figure;
mesh(reshape(x2, Nx, Ny));
title('Numerical Solution - 9-Point Method');

figure;
mesh(reshape(abs(x2 - U_ANAL), Nx, Ny));
title('Absolute Error - 9-Point Method');
t2=toc(t1);
disp(['Συνολικός χρόνος που χρειάστηκε για τους υπολογισμούς:',num2str(t2)]);
