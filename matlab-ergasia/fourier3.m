clear;
clc;

T = 3;  % Period
omega0 = 2 * pi / T;  % Fundamental angular frequency
n_max = 400;  % Number of terms in the Fourier series
t = 0:0.01:T;  % Time vector for plotting

% Function f(t) (piecewise defined)
f = @(t) (t >= 0 & t <= 1) .* t + (t > 1 & t <= 1.5) .* 0 + (t > 1.5 & t <= 2.5) .* (1.5 - t)+(t>2.5 &t<3).*0;

% Fourier coefficients An and Bn
A = @(n) (2/T) * integral(@(t) f(t) .* cos(n * omega0 * t), 0, T);
B = @(n) (2/T) * integral(@(t) f(t) .* sin(n * omega0 * t), 0, T);

% Initialize Fourier series
f_t = zeros(size(t));  % Initialize the Fourier series sum

% Compute Fourier series up to n_max
for n = 1:n_max
    an = A(n);  % Calculate An
    bn = B(n);  % Calculate Bn
    
    % Compute amplitude C_n and phase angle theta_n
    cn = sqrt(an^2 + bn^2);
    theta_n = -atan2(bn, an);  % Phase angle
    
    % Add the nth term to the Fourier series
    f_t = f_t + cn * cos(n * omega0 * t + theta_n);
end

% Plot the result
figure;
plot(t, f_t, 'b-');
xlabel('Time t');
ylabel('f(t)');
title('Reconstruction of f(t) using Fourier Series');
grid on;
% Display a few Fourier coefficients
disp('First few Fourier coefficients (An, Bn):');
for n = 1:5
    fprintf('n = %d, An = %.4f, Bn = %.4f\n', n, A(n), B(n));
    
end
