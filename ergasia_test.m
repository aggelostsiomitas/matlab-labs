clc;
clear;
close all;
n = 4; % Μέγεθος πίνακα

% Ορισμός του πίνακα A
A = [0 0 2 4; 1/2 0 0 0; 0 1/4 0 0; 0 0 1/8 0];
I = eye(n); % Μοναδιαίος πίνακας 4x4

%% (a) Θεώρημα Gerschgorin
gerschgorin_disks = zeros(n, 2); % Πρώτη στήλη κέντρο, δεύτερη ακτίνα

for i = 1:n
    center = A(i,i); % Διαγώνιο στοιχείο
    radius = sum(abs(A(i,:))) - abs(A(i,i)); % Άθροισμα απόλυτων τιμών εκτός διαγωνίου
    gerschgorin_disks(i,:) = [center, radius];
end

disp('Δίσκοι Gerschgorin (Κέντρο, Ακτίνα):');
disp(gerschgorin_disks);

%% (b) Μέθοδος Δύναμης (Power Method)
x = rand(n,1); % Τυχαίο αρχικό διάνυσμα
x = x / norm(x); % Κανονικοποίηση
lambda_old = 0;
tol = 1e-6;
max_iter = 100;

for iter = 1:max_iter
    x_new = A * x;
    lambda_new = norm(x_new, inf);
    x = x_new / lambda_new;
    
    if abs(lambda_new - lambda_old) < tol
        break;
    end
    lambda_old = lambda_new;
end

dominant_eigenvalue = lambda_new;
dominant_eigenvector = x;

disp('Κυρίαρχη ιδιοτιμή (Power Method):');
disp(dominant_eigenvalue);
disp('Αντίστοιχο ιδιοδιάνυσμα:');
disp(dominant_eigenvector);

%% (c) Χειροκίνητος Υπολογισμός Ιδιοτιμών και Ιδιοδιανυσμάτων
% Χαρακτηριστικό Πολυώνυμο: det(A - λI) = 0
% Προσεγγιστικά υπολογίζουμε τις ρίζες του πολυωνύμου.

lambda_vals = linspace(-10, 10, 1000); % Δοκιμάζουμε μια γραμμική περιοχή
char_poly_vals = zeros(size(lambda_vals)); % Αποθήκευση τιμών του χαρακτηριστικού πολυωνύμου

% Υπολογισμός του χαρακτηριστικού πολυωνύμου
for i = 1:length(lambda_vals)
    lambda = lambda_vals(i);
    char_poly_vals(i) = det(A - lambda * I);
end

% Εύρεση σημείων μηδενισμού του πολυωνύμου
zero_crossings = find(diff(sign(char_poly_vals)) ~= 0);
eigenvalues = lambda_vals(zero_crossings);

disp('Ιδιοτιμές του Α (Χειροκίνητος Υπολογισμός):');
disp(eigenvalues);

% Υπολογισμός ιδιοδιανυσμάτων για κάθε ιδιοτιμή
eigenvectors = [];

for i = 1:length(eigenvalues)
    lambda_i = eigenvalues(i);
    % Υπολογισμός ιδιοδιανύσματος για κάθε ιδιοτιμή
    % Λύση του συστήματος (A - λI) x = 0
    [V, D] = eig(A - lambda_i * I);
    
    % Το ιδιοδιάνυσμα θα είναι το μη μηδενικό στοιχείο του συνόλου λύσεων
    % Επιλέγουμε την πρώτη στήλη του πίνακα V ως το ιδιοδιάνυσμα
    eigenvectors = [eigenvectors, V(:,1)];
end

disp('Ιδιοδιανύσματα του Α (Χειροκίνητος Υπολογισμός):');
disp(eigenvectors);

%% (d) Μέθοδος Newton-Raphson για εύρεση ιδιοτιμών
% Ορίζουμε την ανώνυμη συνάρτηση για το χαρακτηριστικό πολυώνυμο και την παράγωγο του
char_poly = @(lambda) det(A - lambda * I);
char_poly_derivative = @(lambda) trace(inv(A - lambda * I));  % Η παράγωγος του χαρακτηριστικού πολυωνύμου

% Αρχικές εκτιμήσεις για τις ιδιοτιμές
lambda_guesses = [0.5, -0.5, 0.1, -0.1]; 
tol = 1e-6;
max_iter = 100;

eigenvalues_NR = zeros(size(lambda_guesses)); 

for i = 1:length(lambda_guesses)
    lambda_old = lambda_guesses(i);
    iter = 0;
    
    while iter < max_iter
        lambda_new = lambda_old - char_poly(lambda_old) / char_poly_derivative(lambda_old);
        
        if abs(lambda_new - lambda_old) < tol
            break;
        end
        lambda_old = lambda_new;
        iter = iter + 1;
    end
    
    eigenvalues_NR(i) = lambda_new;
end

disp('Ιδιοτιμές με Newton-Raphson:');
disp(eigenvalues_NR);

%% Μακροχρόνια συμπεριφορά πληθυσμού
growth_rate = dominant_eigenvalue;
disp('Μακροχρόνια πρόβλεψη για τον πληθυσμό:');
if growth_rate > 1
    disp('Ο πληθυσμός των κανθάρων θα αυξάνεται με την πάροδο του χρόνου.');
elseif growth_rate < 1
    disp('Ο πληθυσμός των κανθάρων θα μειώνεται και θα εξαφανιστεί.');
else
    disp('Ο πληθυσμός των κανθάρων θα παραμείνει σταθερός.');
end
