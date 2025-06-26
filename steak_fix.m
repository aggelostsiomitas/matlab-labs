clear;
clc;

% Parameters
a = 1.35e-7;                % Thermal diffusivity
H = 5.08e-2; L = 20.32e-2;  % Steak dimensions (m)
K = 0.04;                   % Neumann boundary factor
T_initial = 23;             % Initial temperature (°C)
T_pan_bottom = 204;         % Bottom pan temperature (°C)
T_final = 56;               % Desired center temperature (°C)

nx = 60; ny = 60; nt = 3000;
dx = L / nx; dy = H / ny;
dt = 0.2 * min(dx^2, dy^2) / a;   

N = (nx+1)*(ny+1);

% Laplacian Operator
Tx = spdiags(ones(nx+1,1)*[1 -2 1], [-1 0 1], nx+1, nx+1) / dx^2;
Ty = spdiags(ones(ny+1,1)*[1 -2 1], [-1 0 1], ny+1, ny+1) / dy^2;
A = kron(speye(ny+1), Tx) + kron(Ty, speye(nx+1));

I=speye(N);
LHS = I - a * dt * A;

center_x = round((nx+1)/2);
center_y = round((ny+1)/2);

% --- One Flip Simulation ---
fprintf('\n--- Running simulation with 1 flip (Fully Implicit) ---\n');
sim1_start=tic;

T = T_initial * ones(nx+1, ny+1);
T(1,:) = T_pan_bottom;

flip_step = round(nt / 2);
flipped = false;
time_to_final_1flip = NaN;
t1 = tic;

for t = 2:nt
    % Neumann boundaries (top or bottom edge depending on flip)
    cols = 2:nx;
    if ~flipped
        T(end, cols) = T(end-1, cols) + K * dy;
    else
        T(1, cols) = T(2, cols) + K * dy;
    end

    T_vec = reshape(T, [], 1);

    % Apply bottom pan temperature
    if ~flipped
        T_vec(1:nx+1) = T_pan_bottom;
    else
        T_vec(end-nx:end) = T_pan_bottom;
    end

    % Add Neumann left/right boundaries using b
    b = zeros(N,1);
    for j = 1:ny+1
        left_idx  = (j-1)*(nx+1) + 1;
        right_idx = j*(nx+1);
        b(left_idx)  = 2*K/H;
        b(right_idx) = -2*K/H;
    end

    T = reshape(LHS \ (T_vec + b), nx+1, ny+1);

    if t == flip_step && ~flipped
        T = flipud(T);
        flipped = true;
        fprintf('Flipped at t = %.2f s (step %d)\n', t*dt, t);
    end

    T_center = T(center_x, center_y);
    if mod(t, 200) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center, T(1, center_y), T(end, center_y), max(T(:)));
    end

    if T_center >= T_final && isnan(time_to_final_1flip)
        time_to_final_1flip = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center, time_to_final_1flip, t);
        break;
    end
end
sim1_end=toc(sim1_start);

% --- Two Flip Simulation ---
fprintf('\n--- Running simulation with 2 flips (Fully Implicit) ---\n');
sim2_start=tic;

T = T_initial * ones(nx+1, ny+1);
T(1,:) = T_pan_bottom;

first_flip_step = round(nt / 3);
second_flip_step = round(2 * nt / 3);
flipped = false;
second_flip = false;
time_to_final_2flip = NaN;

for t = 2:nt
    % Neumann boundaries (top or bottom edge depending on flip)
    cols = 2:nx;
    if ~flipped
        T(end, cols) = T(end-1, cols) + K * dy;
    else
        T(1, cols) = T(2, cols) + K * dy;
    end

    T_vec = reshape(T, [], 1);

    % Apply bottom pan temperature
    if ~flipped
        T_vec(1:nx+1) = T_pan_bottom;
    else
        T_vec(end-nx:end) = T_pan_bottom;
    end

    % Add Neumann left/right boundaries using b
    b = zeros(N,1);
    for j = 1:ny+1
        left_idx  = (j-1)*(nx+1) + 1;
        right_idx = j*(nx+1);
        b(left_idx)  = 2*K/H;
        b(right_idx) = -2*K/H;
    end

    T = reshape(LHS \ (T_vec + b), nx+1, ny+1);

    if t == first_flip_step && ~flipped
        T = flipud(T);
        flipped = true;
        fprintf('First flip at t = %.2f s (step %d)\n', t*dt, t);
    end

    if t == second_flip_step && flipped && ~second_flip
        T = flipud(T);
        second_flip = true;
        fprintf('Second flip at t = %.2f s (step %d)\n', t*dt, t);
    end

    T_center = T(center_x, center_y);
    if mod(t, 200) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center, T(1, center_y), T(end, center_y), max(T(:)));
    end

    if T_center >= T_final && isnan(time_to_final_2flip)
        time_to_final_2flip = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center, time_to_final_2flip, t);
        break;
    end
end
sim2_end=toc(sim2_start);

% --- Plot Final Results ---
figure;
mesh(T);
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
title(sprintf('Final Temp - 1 Flip (%.2f s)', time_to_final_1flip));
colorbar;

figure;
mesh(T);
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
title(sprintf('Final Temp - 2 Flips (%.2f s)', time_to_final_2flip));
colorbar;

disp(['The difference in the time needed for the 2 simulations is: ', ...
    num2str(abs(sim1_end - sim2_end)), ' seconds']);

if sim1_end < sim2_end
    fprintf('The first simulation was faster by %.2f%%\n', ...
        ((sim2_end - sim1_end) / sim1_end) * 100);
else
    fprintf('The second simulation was faster by %.2f%%\n', ...
        ((sim1_end - sim2_end) / sim2_end) * 100);
end

t2 = toc(t1);
disp(['Total time needed for the entire simulation: ', num2str(t2), ' seconds']);
