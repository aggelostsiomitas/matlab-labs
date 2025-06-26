clear;
clc;

%Definition of differential equations
F=@(x)[4-0.0003*x(1)-0.0004*x(2);2-0.0002*x(1)-0.0001*x(2)];

%Jacobean matrix J(x1,x2)
J=@(x)[-0.0003,-0.0004;-0.0002,-0.0001];

solutions=[0,2000;13333,0];

%Initial guess
x0=[0.5;0.5];

%tolerance
tol=1e-8;

%Maximum number of steps
Nmax=1000;


t1=tic; %starts time
for i=1:Nmax
    %Newton Raphson Method
    x=x0-J(x0)\F(x0);

    if norm(x-x0)<tol &&norm(F(x))<tol
        disp(['Σύγκλιση μετά από ',num2str(i),' επαναλήψεις'])
        break;
    end
    x0=x;
end

if~ismember(x',solutions,'rows')
    solutions=[solutions;x(1),x(2)];
end

%print solutions
disp('Συνολικες λύσεις:');
disp(solutions);

%check if there are other solutions except the two given

if size(solutions,1)>2
    fprintf('\nΥπάρχει και άλλη περίπτωση ισορροπίας.\n')
    disp(solutions(3:end,:));
else
    fprintf('\nΔεν υπάρχουν άλλες περιπτώσεις ισορροπίας.\n');
end

t2=toc(t1); %stop time
disp(['Συνολικός χρόνος που χρειάστηκε:',num2str(t2),'seconds']);

