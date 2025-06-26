clc;
clear;
close all;
n = 4; % Μέγεθος πίνακα
t1=tic;
% Ορισμός του πίνακα A
A = [0 0 2 4; 1/2 0 0 0; 0 1/4 0 0; 0 0 1/8 0];
I = eye(n); 

% (a) Θεώρημα Gerschgorin
gerschgorin_disks = zeros(n, 2); 

for i = 1:n
    center = A(i,i); 
    radius = sum(abs(A(i,:))) - abs(A(i,i)); 
    gerschgorin_disks(i,:) = [center, radius];
end

disp('Δίσκοι Gerschgorin (Κέντρο, Ακτίνα):');
disp(gerschgorin_disks);

% (b) Μέθοδος Δύναμης (Power Method)
x = rand(n,1)+1i*rand(n,1); 
x = x / norm(x); 
lambda_old = 0;
tol = 1e-8;
max_iter = 100;

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

dominant_eigenvalue = lambda_new;
dominant_eigenvector = x;

disp('Κυρίαρχη ιδιοτιμή (Power Method):');
disp(dominant_eigenvalue);
disp('Αντίστοιχο ιδιοδιάνυσμα:');
disp(dominant_eigenvector);

% (c)Υπολογισμός όλων των Ιδιοτιμών και Ιδιοδιανυσμάτων
sigma_vals = [0.5, -0.5, 1i, -1i]; 
eig_vals = zeros(length(sigma_vals), 1); 
eig_vecs = zeros(n, length(sigma_vals)); 

for j = 1:length(sigma_vals)
    sigma = sigma_vals(j); % Shift επιλογή
    v = rand(n,1); 
    v = v / norm(v);
    
    lambda_old = sigma;
    
    for k = 1:max_iter
        w = (A - sigma * eye(n)) \ v; %solve A-σ*Ι=v
        v = w / norm(w);
        lambda = v' * A * v; 
        
        if abs(lambda - lambda_old) < tol
            break;
        end
        lambda_old = lambda;
    end
    eig_vals(j) = lambda; 
    eig_vecs(:, j) = v;   
end

disp('Προσεγγιστικές ιδιοτιμές με τη μέθοδο Shifted Inverse Power:')
disp(eig_vals)

disp('Αντίστοιχα ιδιοδιανύσματα:')
disp(eig_vecs)

% (d) Μέθοδος Newton-Raphson για εύρεση ιδιοτιμών
f = @(x) x.^4 - x/4 - 1/16;
df = @(x) 4*x.^3 - 1/4;

% Set parameters

% Initial guesses for Newton-Raphson
initial_guesses = [ 2, 1+1i,-1-1i,-3-2i];
roots_found = [];

for j = 1:length(initial_guesses)
    x = initial_guesses(j);  
    for iter = 1:max_iter
        fx = f(x);
        dfx = df(x);
        
        % Check if derivative is too small to avoid division errors
        if abs(dfx) < 1e-12
            break;
        end
        
        % Newton-Raphson update
        x_new = x - fx / dfx;
        
        % Convergence check
        if abs(x_new - x) < tol
            break;
        end
        
        x = x_new;
    end
    
    % Check if root is unique before storing
    if all(abs(roots_found - x) > tol)
        roots_found = [roots_found, x];
        
        % Polynomial deflation to find remaining roots
        f = @(x) f(x) / (x - x_new);
        df = @(x) df(x) / (x - x_new);
    end
end

% Display results
fprintf('Roots found:\n');
for k = 1:length(roots_found)
    fprintf('Root %d: %.6f + %.6fi\n', k, real(roots_found(k)), imag(roots_found(k)));
end


% Μακροχρόνια συμπεριφορά πληθυσμού
growth_rate = dominant_eigenvalue;
disp('Μακροχρόνια πρόβλεψη για τον πληθυσμό:');
if growth_rate > 1
    disp('Ο πληθυσμός των κανθάρων θα αυξάνεται με την πάροδο του χρόνου.');
elseif growth_rate < 1
    disp('Ο πληθυσμός των κανθάρων θα μειώνεται και θα εξαφανιστεί.');
else
    disp('Ο πληθυσμός των κανθάρων θα παραμείνει σταθερός.');
end
t2=toc;
disp(t2);
disp(['the time needed needed for all the calculations  was :',num2str(t2),'seconds']);