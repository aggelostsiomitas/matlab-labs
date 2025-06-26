clear;
clc;

T = 3; 
omega0 = 2 * pi / T; 
t = 0:0.0001:T; 
n_max = 200; 
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
    
    cn = sqrt(an^2 + bn^2); 
    theta_n = -atan2(bn, an);
    f_t = f_t + cn * cos(n * omega0 * t + theta_n);
end


Hv = @(t) (1/2) * (dirac(t) - exp(-t / 20) .* heaviside(t));


Hv_vals = (1/2) * (t == 0) - (1/2) * exp(-t / 20) .* (t > 0); 
v0_vals = f_t .* Hv_vals;


omega = linspace(0.1, 100, 1000); 
Hv_freq = (10j .* omega) ./ (1 + 20j .* omega); 


magnitude = abs(Hv_freq);
phase = angle(Hv_freq); 


figure;


subplot(4, 1, 1);
plot(t, f_t, 'b-');
xlabel('Time t');
ylabel('f(t)');
title('Reconstruction of f(t) using Fourier Series');
grid on;


subplot(4, 1, 2);
plot(t, v0_vals, 'k-');
xlabel('Time t');
ylabel('v_0(t)');
title('Ανακατασκευή της v_0(t) = Hν(t) * f(t)');
grid on;


subplot(4, 1, 3);
plot(omega, 20*log10(magnitude), 'b', 'LineWidth', 1.5); 
xlabel('\omega (rad/s)');
ylabel('|H_v(j\omega)| (dB)');
title('Μέτρο της H_v(j\omega)');
grid on;


subplot(4, 1, 4);
plot(omega, phase * (180/pi), 'r', 'LineWidth', 1.5); 
xlabel('\omega (rad/s)');
ylabel('Φάση (°)');
title('Φάση της H_v(j\omega)');
grid on;