clc; clear;

A = [0, 0, 2, 4; 
     1/2, 0, 0, 0; 
     0, 1/4, 0, 0; 
     0, 0, 1/8, 0];

n = size(A,1);
sigma_vals = [0.5, -0.5, 1i, -1i]; 
eig_vals = zeros(length(sigma_vals), 1); 

tol = 1e-6;
max_iter = 100;

for j = 1:length(sigma_vals)
    sigma = sigma_vals(j); % Shift επιλογή
    v = rand(n,1); 
    v = v / norm(v);
    
    lambda_old = sigma;
    
    for k = 1:max_iter
        w = (A - sigma * eye(n)) \ v; % Επίλυση (A - σΙ)x = v
        v = w / norm(w);
        lambda = v' * A * v; % Υπολογισμός ιδιοτιμής
        
        if abs(lambda - lambda_old) < tol
            break;
        end
        lambda_old = lambda;
    end
    eig_vals(j) = lambda; 
end

disp('Προσεγγιστικές ιδιοτιμές με τη μέθοδο Shifted Inverse Power:')
disp(eig_vals)
