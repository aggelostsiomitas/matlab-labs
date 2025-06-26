%%
clear;
clc;
F=@(x) [x(1)^2+x(1)*x(2)-10;...
    x(2)+3*x(1)*x(2)^2-57];

%initial guesss
x0=[0.5;0.5];
G=@(x) [(10-x(1)*x(2))/x(1);...
    (57-x(2))/(3*x(1)*x(2))];

tol=1e-4;
Nmax=10000;
for i=1:Nmax
    x=G(x0);
    %ελεγχω πεδιο τιμων και πεδιο ορισμου γιατι εχω διανυσμα 
    if norm(x-x0)<tol||norm(F(x))<tol
        disp(['convergence after:',num2str(i) 'iterations']);
        break;
    end
    x0=x;
end

%%
clear;
clc;
F=@(x) [x(1)^2+x(1)*x(2)-10;...
    x(2)+3*x(1)*x(2)^2-57];

%initial guesss
x0=[0.5;0.5];
J=@(x) [2*x(1)*x(2),x(1);3*x(2)^2,1+6*x(1)*x(2)];

tol=1e-4;
Nmax=10000;
for i=1:Nmax
    x=x0-J(x0)\F(x0);
    %ελεγχω πεδιο τιμων και πεδιο ορισμου γιατι εχω διανυσμα 
    if norm(x-x0)<tol && norm(F(x))<tol
        disp(['convergence after:',num2str(i) ' iterations']);
        break;
    end
    x0=x;
end
disp(x0);
% ή θα μπορουσα να κανω dx=J(x0)\-F(x0)
% x=x0+dx
% να το φερω δηλαδη σε μορφη 
%J(xn)*δ(xn)=-F(xn), δ(xn)=x(n+!)-x(n)
