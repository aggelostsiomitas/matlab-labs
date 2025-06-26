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
nx = 60; ny = 60;  % Reduced from 100 to 50 for faster computation
dx = L / nx;
dy = H / ny;

% Calculate stable time step
dt = 0.3 * min(dx^2, dy^2) / alpha;  % Adjusted time step for better stability
nt = 1500;  % Increased time steps to ensure enough time for two flips

N = (nx+1)*(ny+1);

% Matrix A with Kronecker products
Ix = speye(nx+1);
Iy = speye(ny+1);
ex = ones(nx+1,1);
ey = ones(ny+1,1);
Tx = spdiags([ex -2*ex ex], -1:1, nx+1, nx+1) / dx^2;
Ty = spdiags([ey -2*ey ey], -1:1, ny+1, ny+1) / dy^2;

A = kron(Iy, Tx) + kron(Ty, Ix);
A_2flip=A;

I=speye(N);
M = I - dt * alpha * A;
M_2flip=M;

% Initialize Temperature matrix for all time steps
T = zeros(N, nt+1);  % Temperature matrix for all time steps
T_2flip = zeros(N, nt+1);  % Temperature matrix for 2-flip simulation

% Set initial conditions
T(:,1) = T_init * ones(N,1);
T_2flip(:,1) = T_init * ones(N,1);

% Apply Dirichlet boundary condition for initial state
T2D = reshape(T(:,1), nx+1, ny+1);
T2D(:,1) = T_oven;
T(:,1) = T2D(:);

T2D_2flip = reshape(T_2flip(:,1), nx+1, ny+1);
T2D_2flip(:,1) = T_oven;
T_2flip(:,1) = T2D_2flip(:);

% Matrix b for Neumann conditions
b = zeros(N,1);
b(nx+1:nx+1:end) = -2 * alpha * K / dx;
b(ny+1:ny+1:end) = 2 * alpha * K / dx;
b(N-(nx+1)+(1:nx+1)) = -2 * alpha * K / dy;

% center of steak 
mid_x = round(nx/2) + 1;
mid_y = round(ny/2) + 1;
mid_index = sub2ind([nx+1, ny+1], mid_x, mid_y);

target_reached = false;
target_reached_2flip = false;

%time for full simulation starts
t1=tic;

t1_1flip = tic; %time for first simulation starts 
fprintf('=====Running simulation for 1 flip====\n');
for i = 2:nt+1
    T(:,i) = M \ (T(:,i-1) + b);
    T2D = reshape(T(:,i), nx+1, ny+1);
    
    % Reapply Dirichlet boundary condition
    T2D(:,1) = T_oven;
    T(:,i) = T2D(:);

    % Get the center temperature
    center_temp = T2D(mid_x, mid_y);

    % Display progress every 20 iterations
    if mod(i, 100) == 0
        disp(['Time step: ', num2str(i), ', Center Temperature: ', num2str(center_temp), '°C']);
    end

    % Stop simulation when center temperature reaches target
    if center_temp >= T_target
        disp(['Center temperature reached the target (', num2str(center_temp), '°C) at time step ', num2str(i)]);
        target_reached = true;
        break;
    end

    % Flip the matrix at halfway point
    if i == round(nt/2)
        disp(['Performing first flip at time step ', num2str(i)]);
        T2D = flipud(T2D);
        T2D(end,:) = T_oven;
        T(:,i) = T2D(:);
    end

    if i == round(2*nt/3)
        disp(['Performing second flip at time step ', num2str(i)]);
        T2D = flipud(T2D);
        T2D(end,:) = T_oven;
        T(:,i) = T2D(:);
    end
end

if ~target_reached
    disp('Warning: Target temperature was not reached within the maximum number of time steps');
end
t1_1flip_end = toc(t1_1flip); %time for first simulation ends
disp(['Time needed for the first simulation: ', num2str(t1_1flip_end), ' seconds']);

%second simulation
fprintf('\n\n=====Running simulation for 2 flips=====\n');
t1_2flip = tic; %time for first simulation starts
for i = 2:nt+1
    T_2flip(:,i) = M_2flip \ (T_2flip(:,i-1) + b);
    T2D_2flip = reshape(T_2flip(:,i), nx+1, ny+1);
    
    % Reapply Dirichlet boundary condition
    T2D_2flip(:,1) = T_oven;
    T_2flip(:,i) = T2D_2flip(:);

    % Get the center temperature
    center_temp = T2D_2flip(mid_x, mid_y);

    % Display progress every 20 iterations
    if mod(i, 100) == 0
        disp(['Time step: ', num2str(i), ', Center Temperature: ', num2str(center_temp), '°C']);
    end

    % Stop simulation when center temperature reaches target
    if center_temp >= T_target
        disp(['Center temperature reached the target (', num2str(center_temp), '°C) at time step ', num2str(i)]);
        target_reached_2flip = true;
        break;
    end

    % First flip at 1/3 of the time
    if i == round(nt/3)
        disp(['Performing first flip at time step ', num2str(i)]);
        T2D_2flip = flipud(T2D_2flip);
        T2D_2flip(end,:) = T_oven;
        T_2flip(:,i) = T2D_2flip(:);
    end

    % Second flip at 2/3 of the time
    if i == round(2*nt/3)
        disp(['Performing second flip at time step ', num2str(i)]);
        T2D_2flip = flipud(T2D_2flip);
        T2D_2flip(end,:) = T_oven;
        T_2flip(:,i) = T2D_2flip(:);
    end
end

if ~target_reached_2flip
    disp('Warning: Target temperature was not reached within the maximum number of time steps');
end
t1_2flip_end = toc(t1_2flip); %time for second simulation ends

disp(['Time needed for the second simulation: ', num2str(t1_2flip_end), ' seconds']);
if t1_2flip_end < t1_1flip_end
    fprintf('The second simulation was faster by %.2f%%\n', ...
        ((t1_1flip_end-t1_2flip_end)/t1_1flip_end)*100);
elseif t1_1flip_end < t1_2flip_end
    fprintf('The first simulation was faster by %.2f%%\n', ...
        ((t1_2flip_end-t1_1flip_end)/t1_2flip_end)*100);
end

% Final Mesh Visualization for first simulation
figure(1);
mesh(T2D);
title('Final Temperature Distribution (1 Flip)');
xlabel('X');
ylabel('Y');
zlabel('Temperature (°C)');

% Final Mesh Visualization for second simulation
figure(2);
mesh(T2D_2flip);
title('Final Temperature Distribution (2 Flips)');
xlabel('X');
ylabel('Y');
zlabel('Temperature (°C)');

t2=toc(t1); %time for full simulation ends
disp(['The total time needed for the full simulation was:',num2str(t2)]);