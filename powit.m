function[y,lambda]=powit(A,Nmax,tol)
N=length(A);
%initial guess
x0=randn(N,1); %βαζυ αυτην την συναρτηση για να παρω τυχαιες τιμες 
% μεσα στο ευρος το οποίο έχω (δεν φευγω εκτος)

%normalize 
x0=x0/norm(x0,inf);

% compute lambda 
lambdax0=(x0'*A*x0)/(x0'*x0);
for i=1:Nmax
    y=A*x0;
    y=y/norm(y,inf);
    lambda=(y'*A*y)/(y'*y);

    if (abs(lambdax0-lambda))<tol
        disp(['Iterations',num2str(i)]);
        break;
    end
    %make the old to new 
    lambdax0=lambda;
    x0=y;
end
y=y/norm(y);
disp(lambda);