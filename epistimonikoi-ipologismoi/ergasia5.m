clear;
clc;

% Parameters
alpha = 1.35e-7;
K = 0.04; % Neumann boundary condition (flux)
T_init = 23; % Initial temperature
T_oven = 204; % Dirichlet boundary temperature
T_target = 56; % Target temperature
L = 0.2032;
H = 0.0508;
nx = 30; ny = 30;
dx = L / nx;
dy = H / ny;
nt = 10000;
dt =alpha/nt;
N = (nx+1)*(ny+1);

% Matrix A with Kronecker products
Ix = speye(nx+1);
Iy = speye(ny+1);
ex = ones(nx+1,1);
ey = ones(ny+1,1);
Tx = spdiags([ex -2*ex ex], -1:1, nx+1, nx+1) / dx^2;
Ty = spdiags([ey -2*ey ey], -1:1, ny+1, ny+1) / dy^2;
A = kron(Iy, Tx) + kron(Ty, Ix);

M = speye(N) - dt * alpha * A;

% Initialize Temperature matrix
T = T_init * ones(N,1);

% Apply Dirichlet boundary condition
T2D = reshape(T, nx+1, ny+1);
T2D(:,1) = T_oven;
T = T2D(:);

% Matrix b for Neumann conditions
b = zeros(N,1);
b(nx+1:nx+1:end) = -2 * alpha * K / dx;
b(ny+1:ny+1:end) = 2 * alpha * K / dx;
b(N-(nx+1)+(1:nx+1)) = -2 * alpha * K / dy;

% Time-Stepping Loop with progress display
mid_x = round(nx/2) + 1;
mid_y = round(ny/2) + 1;
mid_index = sub2ind([nx+1, ny+1], mid_x, mid_y);

for i = 2:nt+1
    T = M \ (T + b);
    T2D = reshape(T, nx+1, ny+1);

    % Get the center temperature
    center_temp = T2D(mid_x, mid_y);

    % Display progress every 100 iterations
    if mod(i, 100) == 0
        disp(['Time step: ', num2str(i), ', Center Temperature: ', num2str(center_temp), '°C']);
    end

    % Stop simulation when center temperature reaches target
    if abs(center_temp - T_target) < 0.1
        disp(['Center temperature reached the target (', num2str(center_temp), '°C) at time step ', num2str(i)]);
        break;
    end

    % Flip the matrix at halfway point
    if i == round(nt/2)
        T2D = flipud(T2D);
        T2D(end,:)=T_oven;
        T = T2D(:);
    end
end

% Final Mesh Visualization
figure;
mesh(T2D);
