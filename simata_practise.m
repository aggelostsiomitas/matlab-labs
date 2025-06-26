% Ορισμός παραμέτρων χρόνου
dt = 0.001; % Βήμα δειγματοληψίας
t_x = -1:dt:1; % Χρονικός άξονας για x(t)
t_h = -1:dt:1; % Χρονικός άξονας για h(t)

% Ορισμός της συνάρτησης x(t)
x = zeros(size(t_x));
x(t_x >= -1 & t_x <= 0) = t_x(t_x >= -1 & t_x <= 0) + 1;
x(t_x > 0 & t_x <= 1) = -t_x(t_x > 0 & t_x <= 1) + 1;

% Ορισμός της συνάρτησης h(t)
h = zeros(size(t_h));
h(t_h >= -1 & t_h <= 1) = t_h(t_h >= -1 & t_h <= 1);

% Εκτέλεση συνέλιξης
y_conv = conv(x, h) * dt; % Πολλαπλασιασμός με dt για προσέγγιση ολοκληρώματος

% Δημιουργία χρονικού άξονα για το αποτέλεσμα
t_y = linspace(-2, 2, length(y_conv));

% Σχεδίαση αποτελεσμάτων
figure;

subplot(3,1,1);
plot(t_x, x, 'LineWidth', 2);
title('Σήμα εισόδου x(t)');
xlabel('Χρόνος (t)');
ylabel('Πλάτος');
xlim([-1.5 1.5]);
grid on;

subplot(3,1,2);
plot(t_h, h, 'LineWidth', 2);
title('Κρουστική απόκριση h(t)');
xlabel('Χρόνος (t)');
ylabel('Πλάτος');
xlim([-1.5 1.5]);
grid on;

subplot(3,1,3);
plot(t_y, y_conv, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
title('Αποτέλεσμα συνέλιξης y(t) = x(t)*h(t)');
xlabel('Χρόνος (t)');
ylabel('Πλάτος');
xlim([-2.5 2.5]);
grid on;

% Βελτιστοποίηση της εμφάνισης
set(gcf, 'Position', [100 100 800 600]);

%%
% Ορισμός της αναλυτικής λύσης
t_analytical = -2:0.01:2;
y_analytical = zeros(size(t_analytical));

for i = 1:length(t_analytical)
    t = t_analytical(i);
    if t >= -2 && t < -1
        y_analytical(i) = (1/6)*(t+2)^3;
    elseif t >= -1 && t < 0
        y_analytical(i) = (-1/2)*t^3 - t^2 + 2/3;
    elseif t >= 0 && t < 1
        y_analytical(i) = (1/2)*t^3 - t^2 + 2/3;
    elseif t >= 1 && t < 2
        y_analytical(i) = (-1/6)*(t-2)^3;
    end
end

% Σύγκριση με αριθμητικό αποτέλεσμα
figure;
plot(t_y, y_conv, 'b', 'LineWidth', 2); hold on;
plot(t_analytical, y_analytical, 'r--', 'LineWidth', 1.5);
legend('Αριθμητική', 'Αναλυτική');
title('Σύγκριση Αποτελεσμάτων');
xlabel('Χρόνος (t)');
ylabel('y(t)');
grid on;

%%
clc; clear; close all;

% Ορισμός σημάτων
dt = 0.001;
t = -3:dt:3;

% Σήμα x(t) (τριγωνικό)
x = zeros(size(t));
x(t >= -1 & t <= 0) = t(t >= -1 & t <= 0) + 1;
x(t > 0 & t <= 1) = -t(t > 0 & t <= 1) + 1;

% Σήμα h(t) (κλίση)
h = zeros(size(t));
h(t >= -1 & t <= 1) = t(t >= -1 & t <= 1);

% Συνέλιξη
y = conv(x, h)*dt;
ty = linspace(-2, 2, length(y));

% Σχεδίαση ανά διαστήματα επικάλυψης
figure('Position', [100, 100, 1200, 800]);

% Αρχικά σήματα
subplot(4,2,[1 2]);
plot(t, x, 'b', 'LineWidth', 2); hold on;
plot(t, h, 'r', 'LineWidth', 2);
title('Αρχικά Σήματα');
legend('x(t)', 'h(t)');
xlim([-2.5 2.5]);
grid on;

% Διάστημα 1: t ∈ [-2, -1) - Μερική επικάλυψη
subplot(4,2,3);
idx = ty >= -2 & ty < -1;
plot(ty(idx), y(idx), 'm', 'LineWidth', 2);
title('Επικάλυψη για t ∈ [-2, -1)');
xlim([-2 -1]);
grid on;

% Διάστημα 2: t ∈ [-1, 0) - Πλήρης επικάλυψη (αριστερά)
subplot(4,2,4);
idx = ty >= -1 & ty < 0;
plot(ty(idx), y(idx), 'g', 'LineWidth', 2);
title('Επικάλυψη για t ∈ [-1, 0)');
xlim([-1 0]);
grid on;

% Διάστημα 3: t ∈ [0, 1) - Πλήρης επικάλυψη (δεξιά)
subplot(4,2,5);
idx = ty >= 0 & ty < 1;
plot(ty(idx), y(idx), 'c', 'LineWidth', 2);
title('Επικάλυψη για t ∈ [0, 1)');
xlim([0 1]);
grid on;

% Διάστημα 4: t ∈ [1, 2) - Μερική επικάλυψη
subplot(4,2,6);
idx = ty >= 1 & ty < 2;
plot(ty(idx), y(idx), 'k', 'LineWidth', 2);
title('Επικάλυψη για t ∈ [1, 2)');
xlim([1 2]);
grid on;

% Ολική συνέλιξη
subplot(4,2,[7 8]);
hold on;
plot(ty(ty >= -2 & ty < -1), y(ty >= -2 & ty < -1), 'm', 'LineWidth', 2);
plot(ty(ty >= -1 & ty < 0), y(ty >= -1 & ty < 0), 'g', 'LineWidth', 2);
plot(ty(ty >= 0 & ty < 1), y(ty >= 0 & ty < 1), 'c', 'LineWidth', 2);
plot(ty(ty >= 1 & ty < 2), y(ty >= 1 & ty < 2), 'k', 'LineWidth', 2);
title('Ολική Συνέλιξη y(t) = x(t)*h(t)');
legend('t ∈ [-2,-1)', 't ∈ [-1,0)', 't ∈ [0,1)', 't ∈ [1,2)');
xlabel('Χρόνος (t)');
ylabel('Πλάτος');
xlim([-2.5 2.5]);
grid on;
hold off;

% Επιπλέον οπτικοποίηση επικαλύψεων
figure('Position', [100, 100, 1200, 600]);
t_shift = -1.5; % Παράδειγμα για t = -1.5

% Δημιουργούμε το σήμα h(t_shift-τ)
h_shifted = zeros(size(t));
shifted_indices = t_shift - t;
valid_indices = shifted_indices >= -1 & shifted_indices <= 1;
h_shifted(valid_indices) = shifted_indices(valid_indices);

% Αρχικά σήματα
subplot(1,2,1);
plot(t, x, 'b', 'LineWidth', 2); hold on;
plot(t, h_shifted, 'r', 'LineWidth', 2);
title(['Σήματα x(τ) και h(',num2str(t_shift),'-τ)']);
xlabel('τ');
legend('x(τ)', ['h(',num2str(t_shift),'-τ)']);
xlim([-3 3]);
grid on;

% Περιοχή επικάλυψης
subplot(1,2,2);
area(t, x.*h_shifted, 'FaceColor', 'y');
title(['Περιοχή Επικάλυψης για t = ',num2str(t_shift)]);
xlabel('τ');
ylabel('x(τ)h(t-τ)');
xlim([-3 3]);
grid on;