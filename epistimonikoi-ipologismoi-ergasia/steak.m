clear; 
clc;

% Parameters
a = 1.35e-7;                % Thermal diffusivity
H = 5.08e-2; L = 20.32e-2;  % Steak dimensions (m)
K = 0.04;                   % Neumann boundary factor
T_initial = 23;            % Initial temperature (°C)
T_pan_bottom = 204;        % Bottom pan temperature (°C)
T_pan_top = 204;           % Top pan temperature after flip (°C)
T_final = 56;              % Desired center temperature (°C)

nx = 60; ny = 60; nt = 2000;
dx = L / nx; dy = H / ny;
dt = 0.2 * min(dx^2, dy^2) / a;

t1=tic;
% 2D Laplacian Operator using Kronecker products
N = (nx+1)*(ny+1);
Tx = spdiags(ones(nx+1,1)*[1 -2 1], [-1 0 1], nx+1, nx+1) / dx^2;
Ty = spdiags(ones(ny+1,1)*[1 -2 1], [-1 0 1], ny+1, ny+1) / dy^2;
A  = kron(speye(ny+1), Tx) + kron(Ty, speye(nx+1));

% Crank-Nicholson matrices
I = speye(N);
Bf = I + (a*dt/2)*A;
Bb = I - (a*dt/2)*A;

% Initial condition
T = T_initial * ones(nx+1, ny+1);
T(1,:) = T_pan_bottom;

T_total = zeros(nx+1, ny+1, nt);
T_total(:,:,1) = T;

flip_step = round(nt / 2);
flipped = false;

center_x = round((nx+1)/2);
center_y = round((ny+1)/2);

% Time-stepping loop
for t = 2:nt
    T_vec = reshape(T, [], 1);

    % Apply Dirichlet BC at top/bottom
    if ~flipped
        T_vec(1:nx+1) = T_pan_bottom;
    else
        T_vec(end-nx:end) = T_pan_top;
    end

    rhs = Bf * T_vec;

    % Apply Neumann BCs (left/right)
for i = 1:ny+1
    base = (i-1)*(nx+1);
    rhs(base+1) = rhs(base+1) + a*dt*K/dx * (T(base+2) - T(base+1));
    rhs(base+nx+1) = rhs(base+nx+1) + a*dt*K/dx * (T(base+nx+1) - T(base+nx));
end

    % Apply Neumann BCs (top or bottom)
    for j = 2:nx
        if ~flipped
            top_idx = ny*(nx+1) + j;
            mid_idx = (ny-1)*(nx+1) + j;
            rhs(top_idx) = rhs(top_idx) + a*dt*K/dy * (T(mid_idx) - T(top_idx));
        else
            bot_idx = j + 1;
            next_idx = bot_idx + nx+1;
            rhs(bot_idx) = rhs(bot_idx) + a*dt*K/dy * (T(next_idx) - T(bot_idx));
        end
    end

    % Solve the system and reshape
    T = reshape(Bb \ rhs, nx+1, ny+1);

    % Flip steak
    if t == flip_step && ~flipped
        T = flipud(T);
        flipped = true;
        fprintf('Flipped at t = %.2f s (step %d)\n', t*dt, t);
    end

    T_total(:,:,t) = T;
    T_center = T(center_x, center_y);

    % Status printout
    if mod(t,100) == 0
        fprintf('Step %d: Center = %.2f°C | Bottom = %.2f°C | Top = %.2f°C | Max = %.2f°C\n', ...
            t, T_center, T(1, center_y), T(end, center_y), max(T(:)));
    end

    if T_center >= T_final
        fprintf('Reached %.1f°C at t = %.2f s (step %d)\n', T_center, t*dt, t);
        break;
    end
end

% Final temperature plot
figure;
mesh(T_total(:,:,t));
xlabel('x'); ylabel('y'); zlabel('Temperature (°C)');
title(sprintf('Final Temperature Distribution (%.2f s)', t*dt));

% Center temperature vs time
T_center_hist = squeeze(T_total(center_x, center_y, 1:t));
time = (1:t) * dt;

figure;
plot(time, T_center_hist, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Center Temp (°C)');
title('Center Temperature vs Time'); grid on;

% Bottom and top surface temperature
T_bot_hist = squeeze(T_total(1, center_y, 1:t));
T_top_hist = squeeze(T_total(end, center_y, 1:t));

figure;
plot(time, T_bot_hist, 'b', 'LineWidth', 2, 'DisplayName', 'Bottom');
hold on;
plot(time, T_top_hist, 'g', 'LineWidth', 2, 'DisplayName', 'Top');
xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Surface Temperatures Over Time');
legend('show'); grid on;

t2=toc(t1);
disp(['Time needed for all the calculations is: ' num2str(t2)]);
