clear;
clc;

% Parameters
a = 1.35e-7; % Thermal diffusivity
H = 5.08e-2; L = 20.32e-2; % Steak dimensions (m)
K = 0.04; % Neumann boundary factor
T_initial = 23; % Initial temperature (°C)
T_pan_bottom = 204; % Bottom pan temperature (°C)
T_pan_top = 204; % Top pan temperature after flip (°C)
T_final = 56; % Desired center temperature (°C)

nx = 60; ny = 60; nt = 3000;
dx = L / nx; dy = H / ny;
dt = 0.2 * min(dx^2, dy^2) / a;

t1 = tic;

% --- Simulation with 1 flip (Neumann on T) ---
fprintf('\n--- Running simulation with 1 flip (Neumann on T) ---\n');
t1_1flip_T = tic;

N = (nx+1)*(ny+1);
Tx = spdiags(ones(nx+1,1)*[1 -2 1], [-1 0 1], nx+1, nx+1) / dx^2;
Ty = spdiags(ones(ny+1,1)*[1 -2 1], [-1 0 1], ny+1, ny+1) / dy^2;
A_1flip_T = kron(speye(ny+1), Tx) + kron(Ty, speye(nx+1));

I = speye(N);
Bf_1flip_T = I + (a*dt/2)*A_1flip_T;
Bb_1flip_T = I - (a*dt/2)*A_1flip_T;

T_1flip_T = T_initial * ones(nx+1, ny+1);
T_1flip_T(1,:) = T_pan_bottom;

T_total_1flip_T = zeros(nx+1, ny+1, nt);
T_total_1flip_T(:,:,1) = T_1flip_T;

flip_step_1flip_T = round(nt / 2);
flipped_1flip_T = false;
center_x = round((nx+1)/2);
center_y = round((ny+1)/2);

time_to_final_1flip_T = NaN;
final_T_1flip_T = [];

% Time-stepping loop (1 flip)
for t = 2:nt
    % Neumann BC left/right
    for j = 1:ny+1
        T_1flip_T(1, j) = T_1flip_T(2, j) + K * dx;
        T_1flip_T(nx+1, j) = T_1flip_T(nx, j) + K * dx;
    end

    % Neumann BC top/bottom
    for i = 2:nx
        if ~flipped_1flip_T
            T_1flip_T(ny+1, i) = T_1flip_T(ny, i) + K * dy;
        else
            T_1flip_T(1, i) = T_1flip_T(2, i) + K * dy;
        end
    end

    T_vec = reshape(T_1flip_T, [], 1);

    % Dirichlet BC
    if ~flipped_1flip_T
        T_vec(1:nx+1) = T_pan_bottom;
    else
        T_vec(end-nx:end) = T_pan_top;
    end

    rhs = Bf_1flip_T * T_vec;
    T_1flip_T = reshape(Bb_1flip_T \ rhs, nx+1, ny+1);

    % Flip
    if t == flip_step_1flip_T && ~flipped_1flip_T
        T_1flip_T = flipud(T_1flip_T);
        flipped_1flip_T = true;
        fprintf('Flipped at t = %.2f s (step %d)\n', t*dt, t);
    end

    T_total_1flip_T(:,:,t) = T_1flip_T;

    T_center_1flip_T = T_1flip_T(center_x, center_y);

    if mod(t, 200) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center_1flip_T, T_1flip_T(1, center_y), T_1flip_T(end, center_y), max(T_1flip_T(:)));
    end

    if T_center_1flip_T >= T_final && isnan(time_to_final_1flip_T)
        time_to_final_1flip_T = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center_1flip_T, time_to_final_1flip_T, t);
        final_T_1flip_T = T_1flip_T;
        break;
    end
end

time_needed_1flip_T = toc(t1_1flip_T);
disp(['Time needed for 1 flip: ', num2str(time_needed_1flip_T)]);

center_temp_hist_1flip_T = squeeze(T_total_1flip_T(center_x, center_y, 1:t));
time_vector_1flip_T = (1:t) * dt;

% --- Simulation with 2 flips (Neumann on T) ---
fprintf('\n--- Running simulation with 2 flips (Neumann on T) ---\n');
t1_2flips_T = tic;

A_2flips_T = A_1flip_T;
Bf_2flips_T = Bf_1flip_T;
Bb_2flips_T = Bb_1flip_T;

T_2flips_T = T_initial * ones(nx+1, ny+1);
T_2flips_T(1,:) = T_pan_bottom;
T_total_2flips_T = zeros(nx+1, ny+1, nt);
T_total_2flips_T(:,:,1) = T_2flips_T;

flip_times_2flips_T = round(linspace(nt / 3, 2 * nt / 3, 2));
flipped_2flips_T = false;
flip_count_2flips_T = 0;

time_to_final_2flips_T = NaN;
final_T_2flips_T = [];

% Time-stepping loop (2 flips)
for t = 2:nt
    for j = 1:ny+1
        T_2flips_T(1, j) = T_2flips_T(2, j) + K * dx;
        T_2flips_T(nx+1, j) = T_2flips_T(nx, j) + K * dx;
    end

    for i = 2:nx
        if flip_count_2flips_T == 0 || flip_count_2flips_T == 2
            T_2flips_T(ny+1, i) = T_2flips_T(ny, i) + K * dy;
        else
            T_2flips_T(1, i) = T_2flips_T(2, i) + K * dy;
        end
    end

    T_vec = reshape(T_2flips_T, [], 1);

    if flip_count_2flips_T == 0 || flip_count_2flips_T == 2
        T_vec(1:nx+1) = T_pan_bottom;
    else
        T_vec(end-nx:end) = T_pan_top;
    end

    rhs = Bf_2flips_T * T_vec;
    T_2flips_T = reshape(Bb_2flips_T \ rhs, nx+1, ny+1);

    if flip_count_2flips_T < 2 && t == flip_times_2flips_T(flip_count_2flips_T + 1)
        T_2flips_T = flipud(T_2flips_T);
        flipped_2flips_T = ~flipped_2flips_T;
        flip_count_2flips_T = flip_count_2flips_T + 1;
        fprintf('Flipped at t = %.2f s (step %d, flip %d/2)\n', t*dt, t, flip_count_2flips_T);
    end

    T_total_2flips_T(:,:,t) = T_2flips_T;

    T_center_2flips_T = T_2flips_T(center_x, center_y);

    if mod(t, 200) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center_2flips_T, T_2flips_T(1, center_y), T_2flips_T(end, center_y), max(T_2flips_T(:)));
    end

    if T_center_2flips_T >= T_final && isnan(time_to_final_2flips_T)
        time_to_final_2flips_T = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center_2flips_T, time_to_final_2flips_T, t);
        final_T_2flips_T = T_2flips_T;
        break;
    end
end

time_needed_2flips_T = toc(t1_2flips_T);
disp(['Time needed for 2 flips: ', num2str(time_needed_2flips_T)]);

center_temp_hist_2flips_T = squeeze(T_total_2flips_T(center_x, center_y, 1:t));
time_vector_2flips_T = (1:t) * dt;

% --- Result Comparison ---
fprintf('\n--- Comparison of Results ---\n');
if ~isnan(time_to_final_1flip_T)
    fprintf('1 Flip: %.2f s to reach %.1f°C\n', time_to_final_1flip_T, T_final);
else
    fprintf('1 Flip did not reach %.1f°C.\n', T_final);
end

if ~isnan(time_to_final_2flips_T)
    fprintf('2 Flips: %.2f s to reach %.1f°C\n', time_to_final_2flips_T, T_final);
else
    fprintf('2 Flips did not reach %.1f°C.\n', T_final);
end

% --- Plot Results ---
figure;
plot(time_vector_1flip_T, center_temp_hist_1flip_T, 'm', 'LineWidth', 2); hold on;
plot(time_vector_2flips_T, center_temp_hist_2flips_T, 'g', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Center Temperature (°C)');
title('Center Temperature vs Time');
legend('1 Flip', '2 Flips');
grid on; hold off;

% --- Mesh Plots ---
if ~isempty(final_T_1flip_T)
    figure;
    [X, Y] = meshgrid(linspace(0, L, nx+1), linspace(0, H, ny+1));
    mesh(X*100, Y*100, final_T_1flip_T);
    xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
    title(sprintf('Final Temperature (1 Flip, %.2f s)', time_to_final_1flip_T));
    colorbar;
end

if ~isempty(final_T_2flips_T)
    figure;
    [X, Y] = meshgrid(linspace(0, L, nx+1), linspace(0, H, ny+1));
    mesh(X*100, Y*100, final_T_2flips_T);
    xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
    title(sprintf('Final Temperature (2 Flips, %.2f s)', time_to_final_2flips_T));
    colorbar;
end

t2 = toc(t1);
disp(['Total simulation time: ', num2str(t2), ' seconds']);


%%
clear;
clc;

% Parameters
a = 1.35e-7; % Thermal diffusivity
H = 5.08e-2; L = 20.32e-2; % Steak dimensions (m)
K = 0.04; % Neumann boundary factor
T_initial = 23; % Initial temperature (°C)
T_pan_bottom = 204; % Bottom pan temperature (°C)
T_pan_top = 204; % Top pan temperature after flip (°C)
T_final = 56; % Desired center temperature (°C)
nx = 60; ny = 60; nt = 1000;
dx = L / nx; dy = H / ny;
dt = 0.3 * min(dx^2, dy^2) / a;t1 = tic;

% --- Simulation with 1 flip (Neumann on T) ---
fprintf('\n--- Running simulation with 1 flip (Neumann on T) ---\n');
t1_1flip_T = tic;
N = (nx+1)*(ny+1);
Tx = spdiags(ones(nx+1,1)*[1 -2 1], [-1 0 1], nx+1, nx+1) / dx^2;
Ty = spdiags(ones(ny+1,1)*[1 -2 1], [-1 0 1], ny+1, ny+1) / dy^2;
A_1flip_T = kron(speye(ny+1), Tx) + kron(Ty, speye(nx+1));
I = speye(N);

% Backward Euler method matrix
B_1flip_T = I - (a*dt)*A_1flip_T;
T_1flip_T = T_initial * ones(nx+1, ny+1);
T_1flip_T(1,:) = T_pan_bottom;
T_total_1flip_T = zeros(nx+1, ny+1, nt);
T_total_1flip_T(:,:,1) = T_1flip_T;
flip_step_1flip_T = round(nt / 2);
flipped_1flip_T = false;
center_x = round((nx+1)/2);
center_y = round((ny+1)/2);
time_to_final_1flip_T = NaN;
final_T_1flip_T = [];

% Time-stepping loop (1 flip - direct Dirichlet BC in system)
for t = 2:nt
    for j = 1:ny+1
        T_1flip_T(1, j) = T_1flip_T(2, j) + K * dx;
        T_1flip_T(nx+1, j) = T_1flip_T(nx, j) + K * dx;
    end

    for i = 2:nx
        if ~flipped_1flip_T
            T_1flip_T(ny+1, i) = T_1flip_T(ny, i) + K * dy;
        else
            T_1flip_T(1, i) = T_1flip_T(2, i) + K * dy;
        end
    end

    T_vec_prev = reshape(T_1flip_T, [], 1);
    rhs = T_vec_prev;
    B_mod = B_1flip_T;

    for i = 1:nx+1
        row_index_bottom = i;
        B_mod(row_index_bottom, :) = 0;
        B_mod(row_index_bottom, row_index_bottom) = 1;
        rhs(row_index_bottom) = T_pan_bottom;
        row_index_top = N - nx + i - 1;

        if flipped_1flip_T
            B_mod(row_index_top, :) = 0;
            B_mod(row_index_top, row_index_top) = 1;
            rhs(row_index_top) = T_pan_top;
        end
    end
    T_vec_next = B_mod \ rhs;
    T_1flip_T = reshape(T_vec_next, nx+1, ny+1);

    if t == flip_step_1flip_T && ~flipped_1flip_T
        T_1flip_T = flipud(T_1flip_T);
        flipped_1flip_T = true;
        fprintf('Flipped at t = %.2f s (step %d)\n', t*dt, t);
    end

    T_total_1flip_T(:,:,t) = T_1flip_T;
    T_center_1flip_T = T_1flip_T(center_x, center_y);

    if mod(t, 100) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center_1flip_T, T_1flip_T(1, center_y), T_1flip_T(end, center_y), max(T_1flip_T(:)));
    end

    if T_center_1flip_T >= T_final && isnan(time_to_final_1flip_T)
        time_to_final_1flip_T = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center_1flip_T, time_to_final_1flip_T, t);
        final_T_1flip_T = T_1flip_T;
        break;
    end
end

time_needed_1flip_T = toc(t1_1flip_T);
disp(['Time needed for 1 flip: ', num2str(time_needed_1flip_T)]);
center_temp_hist_1flip_T = squeeze(T_total_1flip_T(center_x, center_y, 1:t));
time_vector_1flip_T = (1:t) * dt;

% --- Simulation with 2 flips (Neumann on T) ---
fprintf('\n--- Running simulation with 2 flips (Neumann on T) ---\n');
t1_2flips_T = tic;
A_2flips_T = A_1flip_T;
B_2flips_T = I - (a*dt)*A_2flips_T;
T_2flips_T = T_initial * ones(nx+1, ny+1);
T_2flips_T(1,:) = T_pan_bottom;
T_total_2flips_T = zeros(nx+1, ny+1, nt);
T_total_2flips_T(:,:,1) = T_2flips_T;
flip_times_2flips_T = round(linspace(nt / 3, 2 * nt / 3, 2));
flipped_2flips_T = false;
flip_count_2flips_T = 0;
time_to_final_2flips_T = NaN;
final_T_2flips_T = [];

for t = 2:nt
    for j = 1:ny+1
        T_2flips_T(1, j) = T_2flips_T(2, j) + K * dx;
        T_2flips_T(nx+1, j) = T_2flips_T(nx, j) + K * dx;
    end

    for i = 2:nx
        if flip_count_2flips_T == 0 || flip_count_2flips_T == 2
            T_2flips_T(ny+1, i) = T_2flips_T(ny, i) + K * dy;
        else
            T_2flips_T(1, i) = T_2flips_T(2, i) + K * dy;
        end
    end
    T_vec_prev = reshape(T_2flips_T, [], 1);
    rhs = T_vec_prev;
    B_mod = B_2flips_T;

    for i = 1:nx+1
        row_index_bottom = i;
        B_mod(row_index_bottom, :) = 0;
        B_mod(row_index_bottom, row_index_bottom) = 1;
        rhs(row_index_bottom) = T_pan_bottom;
        row_index_top = N - nx + i - 1;

        if flip_count_2flips_T == 1
            B_mod(row_index_top, :) = 0;
            B_mod(row_index_top, row_index_top) = 1;
            rhs(row_index_top) = T_pan_top;
        end
    end

    T_vec_next = B_mod \ rhs;
    T_2flips_T = reshape(T_vec_next, nx+1, ny+1);

    if flip_count_2flips_T < 2 && t == flip_times_2flips_T(flip_count_2flips_T + 1)
        T_2flips_T = flipud(T_2flips_T);
        flipped_2flips_T = ~flipped_2flips_T;
        flip_count_2flips_T = flip_count_2flips_T + 1;
        fprintf('Flipped at t = %.2f s (step %d, flip %d/2)\n', t*dt, t, flip_count_2flips_T);
    end

    T_total_2flips_T(:,:,t) = T_2flips_T;
    T_center_2flips_T = T_2flips_T(center_x, center_y);

    if mod(t, 100) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center_2flips_T, T_2flips_T(1, center_y), T_2flips_T(end, center_y), max(T_2flips_T(:)));
    end

    if T_center_2flips_T >= T_final && isnan(time_to_final_2flips_T)
        time_to_final_2flips_T = t * dt;
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center_2flips_T, time_to_final_2flips_T, t);
        final_T_2flips_T = T_2flips_T;
        break;
    end
end

time_needed_2flips_T = toc(t1_2flips_T);

disp(['Time needed for 2 flips: ', num2str(time_needed_2flips_T)]);
center_temp_hist_2flips_T = squeeze(T_total_2flips_T(center_x, center_y, 1:t));
time_vector_2flips_T = (1:t) * dt;

% --- Result Comparison ---
fprintf('\n--- Comparison of Results ---\n');
if ~isnan(time_to_final_1flip_T)
    fprintf('1 Flip: %.2f s to reach %.1f°C\n', time_to_final_1flip_T, T_final);
else
    fprintf('1 Flip did not reach %.1f°C.\n', T_final);
end
if ~isnan(time_to_final_2flips_T)
    fprintf('2 Flips: %.2f s to reach %.1f°C\n', time_to_final_2flips_T, T_final);
else
    fprintf('2 Flips did not reach %.1f°C.\n', T_final);
end

fprintf('The difference in time between the two simulations is:%.4f\n',abs(time_to_final_2flips_T-time_to_final_1flip_T));
if time_to_final_1flip_T > time_to_final_2flips_T
        fprintf('The second simulation was faster by a percentage of %.4f%%\n', ((time_to_final_1flip_T - time_to_final_2flips_T) / time_to_final_1flip_T) * 100);
    elseif time_to_final_2flips_T > time_to_final_1flip_T
        fprintf('The first simulation was faster by a percentage of %.4f%%\n', ((time_to_final_2flips_T - time_to_final_1flip_T) / time_to_final_2flips_T) * 100);
    else
        fprintf('The two simulations took the same amount of time.\n');
 end
% --- Plot Results ---
figure;
plot(time_vector_1flip_T, center_temp_hist_1flip_T, 'm', 'LineWidth', 2); hold on;
plot(time_vector_2flips_T, center_temp_hist_2flips_T, 'g', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Center Temperature (°C)');
title('Center Temperature vs Time');
legend('1 Flip', '2 Flips');
grid on; hold off;
% --- Mesh Plots ---
if ~isempty(final_T_1flip_T)
    figure;
    [X, Y] = meshgrid(linspace(0, L, nx+1), linspace(0, H, ny+1));
    mesh(X*100, Y*100, final_T_1flip_T);
    xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
    title(sprintf('Final Temperature (1 Flip, %.2f s)', time_to_final_1flip_T));
    colorbar;
end
if ~isempty(final_T_2flips_T)
    figure;
    [X, Y] = meshgrid(linspace(0, L, nx+1), linspace(0, H, ny+1));
    mesh(X*100, Y*100, final_T_2flips_T);
    xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
    title(sprintf('Final Temperature (2 Flips, %.2f s)', time_to_final_2flips_T));
    colorbar;
end
t2 = toc(t1);
disp(['Total simulation time: ', num2str(t2), ' seconds']);
