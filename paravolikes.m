%%
clear;
clc;

%parameters 
c=0.1;
nx=10;
nt=10;
%nt=4200 το οριο για r<1/2
Tmax=0.42;
L=1;
%για c=0.5 θελω nt=42000
%Derived parameters 
h=L/nx;
dt=Tmax/nt; %δt
A=-gallery('tridiag',nx-1)/h^2;
I=speye(nx-1);


%initial condition 
x=linspace(0,L,nx+1)';

%Initialize Temp matrix
T=zeros(nx+1,nt+1);%nx+1 μαζι με τα συνορα , χωρις: nx-1
T(:,1)=sin(pi*x);
%loop event 
for i=2:nt+1
    %T(2:nx,i)=(I+dt*c*A)*T(2:nx,i-1);
    %επειδη εχω τα Boundaries στον Τ εγω δεν τα θελω αρα 
    %εχω 11 σημεια το πρωτο και τελευατιο ειναι τα boundaries 
    %θελω τα υπολοιπα 

    %για emplicit μεθοδο τρεχιε για τα οαντα δνε θελω r<1/2
        T(2:nx,i)=(I-dt*c*A)\T(2:nx,i-1);

end
mesh(T);

%%