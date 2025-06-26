%%
eigenvalues = zeros(n, 1);
eigenvectors = zeros(n, n);
eigenvalues(1) = dominant_lambda;
eigenvectors(:, 1) = dominant_idiodianisma;
    
A_new = A - dominant_lambda * (dominant_idiodianisma * dominant_idiodianisma');


for i = 2:n
    x = rand(n, 1); 
    x = x / norm(x); 
    lambda_old = 0;
    
    for iter = 1:Nmax
        x_new = A_new * x;
        lambda_new = norm(x_new, inf);
        x = x_new / lambda_new; 
        
        if abs(lambda_new - lambda_old) < tol
            break;
        end
        lambda_old = lambda_new;
    end
    
    eigenvalues(i) = lambda_new;
    eigenvectors(:, i) = x;
    
    A_new = A_new - lambda_new * (x * x');
end

disp('Ιδιοτιμές (Μέθοδος Δύναμης για όλες):');
disp(eigenvalues);

disp('Ιδιοδιάνυσμα για κάθε ιδιοτιμή:');
disp(eigenvectors);

%d)
