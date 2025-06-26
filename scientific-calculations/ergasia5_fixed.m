clc; clear;
%Parameters
a = 1.35e-7;  %Diffusion Coefficient 
L = 0.2032;
H = 0.0508;
hx = 0.004;
hy = 0.004;
mx = floor(L/hx)+1; %Points with borders
my = floor(H/hy)+1; %Points with borders
dt = 0.1;
Tmax = 10000;
rx = a * dt / hx^2;
ry = a * dt / hy^2;
K = 0.04;
flip_time1 = 2700;     %Χρόνος πρώτου γυρίσματος για ένα γύρισμα
flip_time1_2 = 1430;   %Χρόνος πρώτου γυρίσματος για δύο γυρίσματα
flip_time2 = 3250;     %Χρόνος δεύτερου γυρίσματος για δύο γυρίσματα
flip_time1_3 = 750;    %Χρόνος πρώτου γυρίσματος για τρία γύρισμα
flip_time2_3 = 2200;   %Χρόνος δεύτερου γυρίσματος για τρία γύρισμα
flip_time3 = 3850;     %Χρόνος τρίτου γυρίσματος για τρία γυρίσματα
T_initial = 23;   %Temperature of The Room
T_pan = 204;  %Temperature of The Pan
T_target = 56;   %Temperature of Medium-Rare

%Number of Reversals
f = 3;

%Matrix of The Temperature of Each Grid Point
T = T_initial * ones(my, mx);

%Central Point Coordinates
x_c = round(mx/2);
y_c = round(my/2);

%Sparse coefficient matrix A(2D Backward Euler Method)
N = mx * my;  %All points of the grid(with borders)
Ix = speye(mx);
Iy = speye(my);
ex = ones(mx,1);
ey = ones(my,1);
Tx = spdiags([rx*ex -2*rx*ex rx*ex], [-1 0 1], mx, mx);
Ty = spdiags([ry*ey -2*ry*ey ry*ey], [-1 0 1], my, my);
Dx = kron(Iy, Tx);
Dy = kron(Ty, Ix);
A = speye(N) - (Dx + Dy);

%Method
for t = 1:Tmax    
    %Boundary Condition Dirichlet
    T(end,:) = T_pan; %Down boundary
    
    %Boundaries Conditions Neumann
    T(:,1)   = T(:,2) - K * hx;     %Left Boundary
    T(:,end) = T(:,end-1) - K * hx; %Right Boundary
    T(1,:)   = T(2,:) - K * hy;     %Up Boundary

    %Solution
    T_vec = reshape(T, N, 1);
    T_new = A \ T_vec;
    T = reshape(T_new, my, mx);
    
    %Implementation of Flips
    if f == 1
        if t == flip_time1
            T = flipud(T);
            disp(['Γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
    elseif f == 2
        if t == flip_time1_2
            T = flipud(T);
            disp(['1ο γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
        if t == flip_time2
            T = flipud(T);
            disp(['2ο γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
    else
        if t == flip_time1_3
            T = flipud(T);
            disp(['1ο γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
        if t == flip_time2_3
            T = flipud(T);
            disp(['2ο γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
        if t == flip_time3
            T = flipud(T);
            disp(['3ο γύρισμα μπριζόλας στα: ' num2str(t*dt) ' δευτερόλεπτα'])
        end
    end
    disp(['Θερμοκρασία κεντρικού σημείου:  ' num2str(T(y_c,x_c)) '  βαθμοί κελσίου']);
    %Temperature Test of Central Point
    if T(y_c, x_c) >= T_target
        disp(['Χρόνος για ψήσιμο medium-rare: ' num2str(t*dt) ' δευτερόλεπτα'])
        break
    end
    %Represantation via 3D Graph For Each Second 
    imagesc(T); colorbar; title(['t = ' num2str(t*dt) 'sec']);
    drawnow;
end