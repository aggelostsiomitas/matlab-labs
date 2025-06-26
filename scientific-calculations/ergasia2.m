clear;
clc;

%size of matrix
n=4;

t1=tic; %start time countdown

%Matrix A
A=[0 0 2 4;1/2 0 0 0; 0 1/4 0 0;0 0 1/8 0];

%matrix with ones in main diagonal
I=speye(n);


%a) Gerschgorin theorem
center=diag(A);
radius=sum(abs(A),2)-abs(center);
gerschgorin_disks=[center,radius];

%results
disp('Δίσκοι Gerschgorin (Κέντρο,Ακτίνα):')
disp(gerschgorin_disks);

%Plot the circles
figure;
hold on;

for i=1:n
    x_center=real(center(i));
    y_center=imag(center(i));
    r=radius(i);

    %create the cricles
    theta=linspace (0,2*pi,100);
    x_circle=x_center+r*cos(theta);
    y_circle=y_center+r*sin(theta);

    plot(x_circle,y_circle); %plot circles
    plot(x_center,y_center,'r.') 
end
title('Gerschgorin Disks');
xlabel('Real Axis');
ylabel('Imaginary Axis');
grid on;
axis equal;
hold off;

%b) Power method

%random vector
x=rand(n,1)+1i*rand(n,1);

%normalize
x=x/norm(x);

lambda_old=0;

%tolerance
tol=1e-8;

%max iterations
max_iter=100;

for iter=1:max_iter
    x_new=A*x;
    norm_x_new=norm(x_new);

    x_new=x_new/norm_x_new; %normalize

    lambda_new=(x_new'*A*x_new)/(x_new'*x_new); %find new lambda
    
    %check for tolerance
    if abs(lambda_new-lambda_old)<tol
        break;
    end
    x=x_new;
    lambda_old=lambda_new;
end

%print the result
dominant_eigenvalue=lambda_new;
dominant_eigenvector=x;

disp('Κυρίαρχη Ιδιοτιμή (Power Method):');
disp(dominant_eigenvalue);
disp('Κυρίαρχο Ιδιοδιάνυσμα(Power Method):');
disp(dominant_eigenvector);


%c)Calculation of all eigenvalues and eignvectors

sigma_vals=[0.5,-0.5,1i,-1i];
eig_vals=zeros(length(sigma_vals),1); %store eigenvalues
eig_vecs=zeros(n,length(sigma_vals));  %store eigenvectors

%Use of shifted method to calculate eigenvalues and eigenvectors
for j=1:length(sigma_vals)
    sigma=sigma_vals(j);
    v=rand(n,1);  %random vector
    v=v/norm(v); %normalize
    
    lambda_old=sigma;
    
    for k=1:max_iter
        %solve the A-σ*Ι=v
        w=(A-sigma*I)\v;

        v=w/norm(w);
        lambda=v'*A*v;
    
        %check for tolerance
        if abs(lambda-lambda_old)<tol
            break;
        end
        lambda_old=lambda;
    end

    %store results
    eig_vals(j)=lambda;
    eig_vecs(:,j)=v;
end

%print results
disp('Ιδιοτιμες με χρήση Shifted Inverse Power Method:');
disp(eig_vals);
disp('Ιδιοδιανύσματα με χρήση Shifted Inverse Power Method:');
disp(eig_vecs);

%d) Newton Raphson method to find eigenvalues and eigenvectors

f=@(x)x.^4-x/4-1/16;
df=@(x)4*x.^3-1/4;

%initial guesses
initial_guesses=[2,1+1i,-1-1i,-3-2i];

%matrix to store solutions
roots_found=[];

for j=1:length(initial_guesses)
    x=initial_guesses(j);
    
    for iter=1:max_iter
        fx=f(x);
        dfx=df(x);

        %check if df is too small
        if(abs(dfx))<1e-12
            break;
        end
        %Newton Raphson
        x_new=x-fx/dfx;

        %Convergence check
        if abs(x_new-x)<tol
            break;
        end
    x=x_new;
    end

    %check if rootis unique
    if all(abs(roots_found-x)>tol)
        roots_found=[roots_found,x];

        %elimination of the factorto find the new root
        f=@(x)f(x)/(x-x_new);
        df=@(x)df(x)/(x-x_new);
    end
end

%print solutions
fprintf('Roots found:\n');
for k=1:length(roots_found)
    fprintf('Root %d:%.6f +%.6fi\n',k,real(roots_found(k)),imag(roots_found(k)));
end

%Result for population
growth_rate=dominant_eigenvalue;
disp('Μακροχρόνια Πρόβλεψη για τον πληθυσμό:')

if growth_rate>1
    disp('Ο πληθυσμός των κανθάρων θα συνεχίσει να αυξάνεται με την πάροδο του χρόνου');
elseif growth_rate<1
    disp('Ο πληθυσμός των κανθάρων θα μειώνεται με την πάροδο του χρόνου');
else
    disp('Ο πληθυσμός των κανθάρων θα παραμείνει σταθερός');
end

t2=toc(t1); %end time
disp(['Συνολικός χρόνος που χρειάστηκε για την προσομοίωση:',num2str(t2)])


   







