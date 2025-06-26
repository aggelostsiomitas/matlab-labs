    % Define the period and fundamental frequency
T = 3;
omega0 = 2 * pi / T;

% Define symbolic variables
syms t n

% Define the piecewise function f(t) symbolically
f = piecewise(...
    (t >= 0) & (t <= 1), t, ...         % f(t) = t for 0 <= t <= 1
    (t > 1) & (t <= 1.5), 0, ...       % f(t) = 0 for 1 < t <= 1.5
    (t > 1.5) & (t <= 2.5), 1.5 - t, ... % f(t) = 1.5 - t for 1.5 < t <= 2.5
    (t > 2.5) & (t < 3), 0, ...        % f(t) = 0 for 2.5 < t < 3
    0);                                % Default value outside 0 <= t < 3

% Define the integrals for A(n) and B(n)
A_integral = (2 / T) * int(f * cos(n * omega0 * t), t, 0, T);
B_integral = (2 / T) * int(f * sin(n * omega0 * t), t, 0, T);

% Display the symbolic results
disp('A(n) = ');
disp(A_integral);

disp('B(n) = ');
disp(B_integral);
