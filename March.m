%%
clear;
clc;

N=5;


A=gallery('tridiag',N);

Neig=5;
tol=1e-6;
Nmax=5000;

lambdas=zeros(Neig,1);
V=zeros(N,Neig);

% copy A
A1=A;
for i=1:Neig
    %power iteration
    [V(:,i),lambdas(i)]=powit(A1,Nmax,tol);
    %weiland iteration
A1=A1-lambdas(i)*V(:,i)*V(:,i)';
    
end




%%
clear;
clc;

N=5;


A=gallery('tridiag',N);

Neig=5;
tol=1e-6;
Nmax=5000;

lambdas=zeros(Neig,1);
V=zeros(N,Neig);

% copy A
A1=A;

for i=1:Neig
    %power iteration
    [V(:,i),lambdas(i)]=powit(A1,Nmax,tol);

    %weiland iteration
A1=A1-lambdas(i)*V(:,i)*V(:,i)';
    
end
%%
clear;
clc;

N=5;


A=gallery('tridiag',N);

Neig=5;
tol=1e-6;
Nmax=5000;

lambdas=zeros(Neig,1);
V=zeros(N,Neig);

% copy A
A1=A;

for i=1:Neig
    %power iteration
    [V(i:end,i),lambdas(i)]=powit(A1,Nmax,tol);
x=A1(1,:)/(lambdas(i)*V(i,i));
    %weiland iteration
A1=A1-lambdas(i)*V(i:end,i)*x;

%αφου πολλαπλασιαζω κα διαρω με το lambdas μπορω απλα να το βγαλω 
% x=A1(1,:)/(*V(i,i));
    %weiland iteration
%A1=A1-V(i:end,i)*x;

%take the remaining A
A1=A1(2:end,2:end);
    
end