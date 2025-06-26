%%
clear;
clc;

%Size of matrix
N=6;

%form the matrix 
A=gallery('tridiag',N);

%tolerance
tol=1e-6;

%maximum allowed iterations
Nmax=500;

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

%%
clear;
clc;

%Size of matrix
N=6;

%form the matrix 
A=gallery('tridiag',N);

%tolerance
tol=1e-6;

%maximum allowed iterations
Nmax=500;

%initial guess
x0=randn(N,1); %βαζυ αυτην την συναρτηση για να παρω τυχαιες τιμες 
% μεσα στο ευρος το οποίο έχω (δεν φευγω εκτος)

%normalize 
x0=x0/norm(x0,inf);

% compute lambda 
lambdax0=(x0'*A*x0)/(x0'*x0);
prelambdax0=lambdax0;
for i=1:Nmax
    y=A*x0;
    y=y/norm(y,inf);
    lambda=(y'*A*y)/(y'*y);

    if (abs((lambdax0-lambda))/prelambdax0)<tol
        disp(['Iterations',num2str(i)]);
        break;
    end
    %make the old to new 
    lambdax0=lambda;
    x0=y;
end
disp(lambda);
%%
clear;
clc;

%Size of matrix
N=6;

%form the matrix 
A=gallery('tridiag',N);

%tolerance
tol=1e-6;

%maximum allowed iterations
Nmax=500;

%initial guess
x0=randn(N,1); %βαζυ αυτην την συναρτηση για να παρω τυχαιες τιμες 
% μεσα στο ευρος το οποίο έχω (δεν φευγω εκτος)

%normalize 
x0=x0/norm(x0,inf);

% compute lambda 
lambdax0=(x0'*A*x0)/(x0'*x0);
lambda=lambdax0;
for i=1:Nmax
    y=A\x0;
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
disp(lambda);
%%
clear;
clc;

%Size of matrix
N=6;

%form the matrix 
A=gallery('tridiag',N);

%tolerance
tol=1e-6;

%eigenvalue guess
%sigma=2;

% identity
I=speye(N);

%maximum allowed iterations
Nmax=5000;

%initial guess
x0=randn(N,1); %βαζυ αυτην την συναρτηση για να παρω τυχαιες τιμες 
% μεσα στο ευρος το οποίο έχω (δεν φευγω εκτος)

%normalize 
x0=x0/norm(x0,inf);

% compute lambda 
lambdax0=(x0'*A*x0)/(x0'*x0);
prelmbdax0=lambdax0;
sigma=lambdax0;
for i=1:Nmax
    y=(A-sigma*I)\x0;
    y=y/norm(y,inf);
    lambda=(y'*A*y)/(y'*y);
    sigma=lambda;
    if (abs(lambdax0-lambda))<tol
        disp(['Iterations ',num2str(i)]);
        break;
    end
    %make the old to new 
    lambdax0=lambda;
    x0=y;
end
disp(lambda);
