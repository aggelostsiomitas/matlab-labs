%%
clc;
clear;

%parameters
c=0.5;
nx=100;
nt=24000;
th=0.5;
Tmax=5*0.42;
dt=Tmax/nt;
L=1;
h=L/nx;

A=-gallery('tridiag',nx)/h^2;
I=speye(nx);
A(1,2)=2/h^2;
x=linspace(0,L,nx+1)';

T=zeros(nx+1,nt+1);%nx+1 μαζι με τα συνορα , χωρις: nx-1
T(:,1)=cos(pi*x/2);
%loop event 
for i=2:nt+1
    T(1:nx,i)=(I+dt*c*A)*T(1:nx,i-1);
    %επειδη εχω τα Boundaries στον Τ εγω δεν τα θελω αρα 
    %εχω 11 σημεια το πρωτο και τελευταιο ειναι τα boundaries 
    %θελω τα υπολοιπα 

    %για emplicit μεθοδο τρεχει για οταν δεν θελω r<1/2
       % T(2:nx,i)=(I-dt*c*A)\T(2:nx,i-1);

end
mesh(T);

%%
%με δυο neymann
clc;
clear;

%parameters
c=0.5;
nx=100;
nt=24000;
th=0.5;
Tmax=5*0.42;
dt=Tmax/nt;
L=1;
h=L/nx;

A=-gallery('tridiag',nx+1)/h^2;
I=speye(nx+1);
A(1,2)=2/h^2;
A(nx+1,nx)=2/h^2;
x=linspace(0,L,nx+1)';

T=zeros(nx+1,nt+1);%nx+1 μαζι με τα συνορα , χωρις: nx-1
T(:,1)=cos(pi*x);
%loop event 
for i=2:nt+1
    T(1:nx+1,i)=(I+dt*c*A)*T(1:nx+1,i-1);
  
end
mesh(T);
%%
clc;
clear;

%parameters
c=0.5;
nx=100;
nt=24000;
th=0.5;
Tmax=5*0.42;
dt=Tmax/nt;
L=1;
h=L/nx;
a=1;

A=-gallery('tridiag',nx+1)/h^2;
I=speye(nx+1);
A(1,2)=2/h^2;
A(nx+1,nx)=2/h^2;
x=linspace(0,L,nx+1)';

T=zeros(nx+1,nt+1);%nx+1 μαζι με τα συνορα , χωρις: nx-1
T(:,1)=cos(pi*x);

b=zeros(nx+1,1);
b(nx+1)=2/h*a;

%loop event 
for i=2:nt+1
    T(1:nx+1,i)=(I+dt*c*A)*T(1:nx+1,i-1)+b;
  
end
mesh(T);

%%
clc;
clear;

%parameters
c=0.5;
nx=100;
nt=24000;
th=0.5;
Tmax=5*0.42;
dt=Tmax/nt;
L=1;
h=L/nx;
a=-1;

A=-gallery('tridiag',nx+1)/h^2;
I=speye(nx+1);
A(1,2)=2/h^2;
A(nx+1,nx)=2/h^2;
x=linspace(0,L,nx+1)';

T=zeros(nx+1,nt+1);%nx+1 μαζι με τα συνορα , χωρις: nx-1
T(:,1)=cos(pi*x);

b=zeros(nx+1,1);
b(nx+1)=2/h*a;

%loop event 
for i=2:nt+1
    T(1:nx+1,i)=(I-dt*c*A)\(T(1:nx+1,i-1)+ b);
  
end
mesh(T);