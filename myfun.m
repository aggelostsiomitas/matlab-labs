function [lambdas, v_init] = myfun(A, N, x, v_init)
    % Υπολογισμός ιδιοτιμών με τη μέθοδο της δύναμης για τριδιαγώνιο πίνακα
    % A: τριδιαγώνιος πίνακας
    % N: διάσταση του πίνακα
    % x: αρχικό διάνυσμα
    % v_init: αρχικό ιδιοδιάνυσμα
    
    num_eigenvalues = N; % Αριθμός ιδιοτιμών προς υπολογισμό
    lambdas = zeros(num_eigenvalues, 1); % Αποθήκευση ιδιοτιμών

    tol = 1e-10; % Όριο σύγκλισης
    max_iters = 1000; % Μέγιστος αριθμός επαναλήψεων

    % Αρχικοποίηση του A1 ως το A
    A1 = A;

    for i = 1:num_eigenvalues
        x = x / norm(x); % Κανονικοποίηση
        lambda_old = 0;
        
        % Υπολογισμός A*x και αποθήκευση στο y
        y = A1 * x;

        for k = 1:max_iters
            lambda_new = x' * y; % Υπολογισμός ιδιοτιμής
            x = y / norm(y); % Κανονικοποίηση

            % Έλεγχος σύγκλισης
            if abs(lambda_new - lambda_old) < tol
                break;
            end
            lambda_old = lambda_new;

            % Ενημέρωση του y χωρίς να επαναϋπολογίζουμε το A
            y = A1 * x;
        end

        % Αποθήκευση της ιδιοτιμής
        lambdas(i) = lambda_new;
        
        
        % Ενημέρωση του v_init (το ιδιοδιάστημα)
        v_init = x;
        
        disp(['Iterations',num2str(i)]);
    disp(lambda_new);
        % Ενημέρωση του πίνακα A1 με την εξίσωση Wieland iteration
        A1 = A1 - lambda_new * (x * (v_init' * x)); % WIELAND iteration
    end
    
end
