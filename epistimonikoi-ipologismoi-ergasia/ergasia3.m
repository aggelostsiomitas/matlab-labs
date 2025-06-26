clc;
clear;
t1=tic;
A = [0.1, 0.2, 0.2, 0.6, 0.2;
     0.1, 0.1, 0.1, 0.1, 0.2;
     0.1, 0.3, 0.4, 0.1, 0.2;
     0.3, 0.3, 0.1, 0.1, 0.2;
     0.4, 0.1, 0.2, 0.1, 0.2];

n = size(A,1);
x = rand(n,1); 
x = x / norm(x);

tol = 1e-10;
max_iter = 1000;
lambda_old = 0;

for iter = 1:max_iter
    x_new = A * x;
    x_new = x_new / norm(x_new);
    
    lambda_new = (x_new' * A * x_new) / (x_new' * x_new);

    if abs(lambda_new - lambda_old) < tol
        break;
    end
    
    x = x_new;
    lambda_old = lambda_new;
end

pi_power = x / sum(x);

disp('Σταθερή κατανομή (Μέθοδος Δύναμης):');
disp(pi_power');

% Υπολογισμός κέρδους
total_sales = 60e6;
profit_initial = pi_power(1) * total_sales * 12;

fprintf('Ετήσιο κέρδος πριν: %.2f εκατ. ευρώ\n', profit_initial / 1e6);

% Τροποποιημένος πίνακας με νέα πρώτη στήλη
A_mod = A;
A_mod(:,1) = [0.3; 0.1; 0.1; 0.2; 0.3];

x = rand(n,1);
x = x / norm(x);
lambda_old = 0;

for iter = 1:max_iter
    x_new = A_mod * x;
    x_new = x_new / norm(x_new);
    
    lambda_new = (x_new' * A_mod * x_new) / (x_new' * x_new);

    if abs(lambda_new - lambda_old) < tol
        break;
    end
    
    x = x_new;
    lambda_old = lambda_new;
end

pi_mod_power = x / sum(x);

disp('Σταθερή κατανομή μετά την αλλαγή (Μέθοδος Δύναμης):');
disp(pi_mod_power');

profit_modified = pi_mod_power(1) * total_sales * 12;

fprintf('Ετήσιο κέρδος μετά: %.2f εκατ. ευρώ\n', profit_modified / 1e6);
fprintf('Κόστος διαφημιστικής υπηρεσίας: 40 εκατ. ευρώ\n');

profit_gain = profit_modified - profit_initial;
fprintf('Πρόσθετο ετήσιο κέρδος: %.2f εκατ. ευρώ\n', profit_gain / 1e6);

if profit_gain> 40e6
    fprintf('Η επένδυση αξίζει!\n');
else
    fprintf('Η επένδυση δεν αξίζει!\n');
end

t2=toc(t1);
disp(['Συνολικός χρόνος που χρειάστηκε για τους υπολογισμούς:', num2str(t2)]);