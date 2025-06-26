clear;
clc;
format long;

a = 0;
b = 1;
n = 21; % Ensure that n is a multiple of 3
h = (b - a) / n;

f = @(x) 1 / (1 + 2*x + 3*x.^2); % Use anonymous function

tic;
Result = f(a) + f(b);

for i = 1:n-1
    x = a + i * h; % Correct update of x
    
    if mod(i, 3) == 0
        Result = Result + 2 * f(x);
    else
        Result = Result + 3 * f(x);
    end
end

I = ((3 * h) / 8) * Result;
toc;

disp(I);
