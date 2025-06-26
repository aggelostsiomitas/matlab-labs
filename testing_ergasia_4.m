% Ορισμός παραμέτρων
Nx = 9; % Αριθμός σημείων στο x (0.4 / 0.05) +1
Ny = 9; % Αριθμός σημείων στο y
hx = 0.05;
hy = 0.05;
pi_val = pi;
A = zeros(Nx*Ny, Nx*Ny);
b = zeros(Nx*Ny, 1);

% Δόμηση πλέγματος και πίνακα
for i = 2:Nx-1
    for j = 2:Ny-1
        k = (i-1)*Ny + j;
        A(k, k) = -4 / (hx * hy);
        A(k, k-1) = 1 / (hx * hy);
        A(k, k+1) = 1 / (hx * hy);
        A(k, k-Ny) = 1 / (hx * hy);
        A(k, k+Ny) = 1 / (hx * hy);
        
        % Δεξιά πλευρά της εξίσωσης
        x = (i-1) * hx;
        y = (j-1) * hy;
        b(k) = -25 * sin(5 * pi_val / 2 * x) * sin(5 * pi_val / 2 * y);
    end
end

% Εφαρμογή οριακών συνθηκών
for i = 1:Nx
    for j = [1, Ny]
        k = (i-1)*Ny + j;
        A(k, :) = 0;
        A(k, k) = 1;
        b(k) = 0;
    end
end

for j = 1:Ny
    for i = [1, Nx]
        k = (i-1)*Ny + j;
        A(k, :) = 0;
        A(k, k) = 1;
        b(k) = 0;
    end
end

% Επίλυση του συστήματος
u = A \ b;

% Αναμόρφωση αποτελεσμάτων και απεικόνιση
% Αναμόρφωση αποτελεσμάτων και απεικόνιση με mesh
U_matrix = reshape(u, [Ny, Nx]);
[X, Y] = meshgrid(0:hx:0.4, 0:hy:0.4);
mesh(X, Y, U_matrix);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Αριθμητική επίλυση της μερικής διαφορικής εξίσωσης');
colorbar;


