clear;
clc;
a = [3 2 -1 2 -2; 1 4 0 2 1; 2 1 2 -1 3; 1 1 -1 3 4];
n = 4;

for k = 1:n-1               % Loop over each pivot column
    for i = k+1:n           % Loop over the rows below the pivot
        m = a(i, k) / a(k, k);  % Compute the multiplier (store in variable m)
        for j = k:n+1       % Update the remaining elements in the row
            a(i, j) = a(i, j) - m * a(k, j);
        end      % Explicitly set the lower triangular element to zero
    end
end

disp(a); % Display the resulting matrix
