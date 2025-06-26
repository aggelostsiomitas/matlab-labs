%%
clear;
clc;

% Parameters
a = 1.35e-7;  % Thermal diffusivity
H = 5.08e-2; L = 20.32e-2;  % Steak dimensions (m)
K = 0.04;  % Neumann boundary factor
T_initial = 23;  % Initial temperature (°C)
T_pan_bottom = 204;  % Fixed bottom temperature (°C)
T_pan_top = 204;  % Top pan temperature after flip (°C)
T_final = 56;  % Desired center temperature (°C)

% Grid parameters
nx = 40; ny = 40; nt = 2000;
hx = L / nx; hy = H / ny;
dt = 0.1 * min(hx^2, hy^2) / a;  % Stable time step

% Create tridiagonal matrices for spatial derivatives
Ax = spdiags([-ones(nx+1,1), 2*ones(nx+1,1), -ones(nx+1,1)], -1:1, nx+1, nx+1) / hx^2;
Ay = spdiags([-ones(ny+1,1), 2*ones(ny+1,1), -ones(ny+1,1)], -1:1, ny+1, ny+1) / hy^2;

% Identity matrices
Ix = speye(nx+1);
Iy = speye(ny+1);

% Construct sparse Laplacian for 2D problem
A = kron(Iy, Ax) + kron(Ay, Ix);

% Initialize solution matrix
T = T_initial * ones(nx+1, ny+1);

% Initialize boundary condition vector
b = zeros((nx+1)*(ny+1), 1);

% Apply Neumann boundary conditions
% Left boundary
b(1:ny+1) = 2*K/hx*a;

% Right boundary
b(end-ny:end) = -2*K/hx*a;

% Top boundary
b(1:nx+1:end) = 2*K/hy*a;

% Precompute matrix for solving
M = speye(size(A)) - dt*a*A;

% Time-stepping loop
for i = 2:nt+1
    % Flip the steak at the halfway point
    if i == round(nt/2)
        T = flipud(T);
        % Update boundary conditions after flip
        T(end,:) = T_pan_bottom;  % Bottom becomes new top
        T(1,:) = T_pan_top;      % Top becomes new bottom
    end
    
    % Apply Dirichlet boundary conditions
    T(end,:) = T_pan_bottom;  % Bottom
    T(1,:) = T_pan_top;      % Top
    
    % Reshape T to a column vector
    T_vec = T(:);
    
    % Solve for next time step
    T_vec = M \ (T_vec + dt*b);
    
    % Reshape back to 2D
    T = reshape(T_vec, nx+1, ny+1);
    
    % Check center temperature
    center_temp = T(round(ny/2), round(nx/2));
    
    % Display progress
    if mod(i, 100) == 0
        disp(['Time step: ', num2str(i), ', Center Temperature: ', num2str(center_temp), '°C']);
    end
    
    % Check if target temperature is reached
    if abs(center_temp - T_final) < 0.1
        disp(['Center temperature reached target (', num2str(center_temp), '°C) at time step ', num2str(i)]);
        break;
    end
end

% Visualize final temperature distribution
figure;
surf(T);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (°C)');
title('Temperature Distribution at Final Time Step');
colorbar;
