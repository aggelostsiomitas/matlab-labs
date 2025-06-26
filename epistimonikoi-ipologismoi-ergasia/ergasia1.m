clear;
clc;
% Ορισμός των διαφορικών εξισώσεων
F = @(x) [4 - 0.0003*x(1) - 0.0004*x(2);
          2 - 0.0002*x(1) - 0.0001*x(2)];
% Ιακωβιανός πίνακας J(x1, x2)
J = @(x) [ -0.0003,-0.0004;-0.0002,-0.0001];
solutions = [0, 2000; 13333, 0];
x0 = [0.5; 0.5];
tol = 1e-4;
Nmax = 10000;

t1=tic; % Start timing before the loop
for i = 1:Nmax
    x = x0 - J(x0) \ F(x0);

    if norm(x - x0) < tol && norm(F(x)) < tol
        disp(['Σύγκλιση μετά από ', num2str(i), ' επαναλήψεις']);
        break;
    end
    x0 = x;
end
if ~ismember(x', solutions, 'rows')
    solutions = [solutions; x(1), x(2)];
end
disp('Συνολικές λύσεις:');
disp(solutions);
if size(solutions, 1) > 2
    fprintf('\nΥπάρχει και άλλη περίπτωση ισορροπίας.\n');
    disp(solutions(3:end, :));
else
    fprintf('\nΔεν υπάρχουν άλλες περιπτώσεις ισορροπίας.\n');
end
t2=toc(t1); % Stop timing after the loop
disp(['Total time needed for the simulation:',num2str(t2),' seconds']);