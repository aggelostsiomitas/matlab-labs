clear;
clc;

T = 3; 
omega0 = 2 * pi / T; 
t = 0:0.0001:T; 
n_max = 200; % Number of Fourier terms
%στους 200 όρους έχει σχηματίζεται η f(t) με μικρές διακυμάνσεις στα σημεία
%αλλαγής , από εδώ και πέρα η αύξηση των όρων θα κάνει αυτές τις
%διακυμάνσεις όλο και λιγότερο διακριτές 

A = @(n) (1/(pi*n))*(+sin(2*pi*n/3)-sin(5*pi*n/3)) + ...
         (3/(2*pi^2*n.^2))*(cos(2*pi*n/3)-1+cos(pi*n)-cos(5*pi*n/3));
B = @(n) (1/(pi*n))*(-cos(2*pi*n/3)+cos(5*pi*n/3)) + ...
         (3/(2*pi^2*n^2))*(+sin(2*pi*n/3)-sin(5*pi*n/3));

f_t = zeros(size(t));
for n = 1:n_max
    an = A(n);
    bn = B(n);
    Cn = sqrt(an^2 + bn^2);
    theta_n = -atan2(bn, an);
    f_t = f_t + Cn * cos(n * omega0 * t + theta_n);
end

H_freq = @(omega) (10j * omega) ./ (1 + 20j * omega);
omega = linspace(0.1, 100, 1000); 
H_vals = H_freq(omega);
magnitude = abs(H_vals);
phase = angle(H_vals);

v0_t = zeros(size(t));
for n = 1:n_max
    an = A(n);
    bn = B(n);
    Cn = sqrt(an^2 + bn^2);
    theta_n = -atan2(bn, an);
    
    omega_n = n * omega0;
    H_val = H_freq(omega_n);
    H_mag = abs(H_val);
    H_phase = angle(H_val);
    
    Cn_prime = H_mag * Cn;
    theta_n_prime = H_phase + theta_n;
    
    v0_t = v0_t + Cn_prime * cos(n * omega0 * t + theta_n_prime);
end

figure;

subplot(3, 1, 1);
plot(t, f_t, 'b-');
xlabel('Time (t)');
ylabel('f(t)');
title('Reconstruction of f(t) using Fourier Series');
grid on;

subplot(3, 1, 2);
yyaxis left;
plot(omega, 20*log10(magnitude), 'b', 'LineWidth', 1.5);
ylabel('|H_v(j\omega)| (dB)');
xlabel('\omega (rad/s)');
title('ΜΕΤΡΟ ΚΑΙ  ΦΑΣΗ ΤΗΣ H_v(j\omega)');
grid on;

yyaxis right;
plot(omega, phase * (180/pi), 'r', 'LineWidth', 1.5);
ylabel('ΦΑΣΗ (°)');

subplot(3, 1, 3);
plot(t, v0_t, 'r-');
xlabel('Time (t)');
ylabel('v_0(t)');
title('Reconstruction of v_0(t) using Transformed Fourier Series');
grid on;
