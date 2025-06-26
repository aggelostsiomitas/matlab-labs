clear;
clc;

t1=tic; %start time

%matrix A
A=[0.1,0.2,0.2,0.6,0.2;
   0.1,0.1,0.1,0.1,0.2;
   0.1,0.3,0.4,0.1,0.2;
   0.3,0.3,0.1,0.1,0.2;
   0.4,0.1,0.2,0.1,0.2];

%copy of A
A_mod=A;

%size A
n=size(A,1);

%random vector
x=rand(n,1);
%normalize 
x=x/norm(x);

%tolerance
tol=1e-10;

%maximum number of iterations 
max_iter=1000;

%initial lambda
lambda_old=0;

%Use of Power method 
for iter=1:max_iter
    x_new=A*x;
    x_new=x_new/norm(x_new); %normalize

    lambda_new=(x_new'*A*x_new)/(x_new'*x_new);

    if abs(lambda_new-lambda_old)<tol
        break;
    end
    x=x_new;
    lambda_old=lambda_new;
end

pi_power=x/sum(x); %normalize the eigenvector
disp('Σταθερή κατανομή(Μέθοδος Δύναμης):');
disp(pi_power');

%profit calculation
total_sales=60e6; %total monthly sales
profit_initial=pi_power(1)*total_sales*12;

fprintf('Ετήσιο κέρδος πριν:%.2f εκατ. ευρώ\n',profit_initial/1e6);

%New matrix A with investment using A coppy 

%use of copy of A
A_mod(:,1)=[0.3,0.1,0.1,0.2,0.3];

%repeating Power method to find the new eigenvector
x=rand(n,1);
x=x/norm(x);
lambda_old=0;
for iter=1:max_iter
     x_new=A_mod*x;
    x_new=x_new/norm(x_new); %normalize

    lambda_new=(x_new'*A_mod*x_new)/(x_new'*x_new);

    if abs(lambda_new-lambda_old)<tol
        break;
    end
    x=x_new;
    lambda_old=lambda_new;
end

pi_mod_power=x/sum(x);

disp('Σταθερή κατανομή μετά την επένδυση(Μέθοδος Δύναμης):')
disp(pi_mod_power');

%calculate the new profit 
profit_modified=pi_mod_power(1)*total_sales*12;

fprintf('Ετήσιο κέρδος μετά: %.2f εκατ. ευρώ\n',profit_modified/1e6);
fprintf('Κόστος διαφήμισηςυπηρεσίας: 40 εκατ. ευρώ\n');

%check if investment is bad or not
profit_gain=profit_modified-profit_initial;
fprintf('Πρόσθετο ετήσιο κέρδος: %.2f εκατ. ευρώ\n',profit_gain/1e6);
if profit_gain>40e6
    fprintf('Η επένδυση αξίζει\n');
else 
    fprintf("Η επένυση δεν αξίζει\n");
end

t2=toc(t1); %end time
disp(['Συνολικος χρόνος που χρειάστηκε:',num2str(t2)]);



















