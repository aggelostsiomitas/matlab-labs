%%
clear;
clc;

%number of intervals 
n=300;
nt=400;
c=0.26;

%mesh size
h=1/n;

%Laplacian
A=-gallery('tridiag',n-1); 
I=speye(n-1);

%time step 
Tmax=4;
dt=Tmax/nt;


%initial condition 
x=linspace(0,1,n+1)';

u=zeros(n+1,n+1);
u(:,1)=sin(3*pi*x);

% compute r
r=c*dt/h;

%first step for t=0
u(2:n,2)=(2*I-r^2*A)\(2*u(2:n,1));

%the remaining steps 
for j=3:nt+1
    u(2:n,j)=(I-r^2*A)\(2*u(2:n,j-1)-u(2:n,j-2));
    plot(x,u(:,j));
    axis([0,1,-1,1]);
    MM=getframe;
end

%plot 
mesh(u);
%εδω βλεπω οτι αποσβαινει το οποιο δεν θα επρεπε να συμβει
%%
%πως θα επρεπε αν ειναι 
clear;
clc;

%number of intervals 
n=300;
nt=400;
c=0.26;

%mesh size
h=1/n;

%Laplacian
A=-gallery('tridiag',n-1); 
I=speye(n-1);

%time step 
Tmax=4;
dt=Tmax/nt;


%initial condition 
x=linspace(0,1,n+1)';

u=zeros(n+1,n+1);
u(:,1)=sin(3*pi*x);

% compute r
r=c*dt/h;

%first step for t=0
u(2:n,2)=(2*I-0.5*r^2*A)\((2*I+0.5*r^2*A)*u(2:n,1));

%the remaining steps 
for j=3:nt+1
    u(2:n,j)=(I-0.25*r^2*A)\((0.5*r^2*A+2*I)*u(2:n,j-1)+(0.25*r^2*A-I)*u(2:n,j-2));
    plot(x,u(:,j));
    axis([0,1,-1,1]);
    MM=getframe;
end

%plot 
mesh(u);

%%
clear;
clc;

%number of intervals 
n=200;
nt=400;
c=0.13;

%mesh size
h=1/n;

%Laplacian
A=-gallery('tridiag',n); 

A(end-1,end)=2;

I=speye(n);

%time step 
Tmax=4;
dt=Tmax/nt;


%initial condition 
x=linspace(0,1,n+1)';

u=zeros(n+1,n+1);
sigma=0.01;
u(:,1)=exp(-(x-0.5).^2/(2*sigma^2)); %dchrislet

% compute r
r=c*dt/h;

%first step for t=0

u(2:n+1,2)=0.5*(2*I+r^2*A)*u(2:n+1,1);

%the remaining steps 
for j=3:nt+1
    u(2:n+1,j)=(2*I+r^2*A)*u(2:n+1,j-1)-u(2:n+1,j-2);
    plot(x,u(:,j));
    axis([0,1,-1,1]);
    MM=getframe;
end

%plot 
mesh(u);

%%
clear;
clc;

%number of intervals 
n=200;
nt=400;
c=0.13;

%mesh size
h=1/n;

%Laplacian
A=-gallery('tridiag',n); 

A(end-1,end)=2;

I=speye(n);

%time step 
Tmax=4;
dt=Tmax/nt;


%initial condition 
x=linspace(0,1,n+1)';

u=zeros(n+1,n+1);
sigma=0.01;
u(:,1)=exp(-(x-0.5).^2/(2*sigma^2)); %dchrislet

% compute r
r=c*dt/h;

%first step for t=0

u(2:n+1,2)=0.5*(2*I+r^2*A)*u(2:n+1,1);

%the remaining steps 
for j=3:nt+1
    u(2:n+1,j)=(2*I+r^2*A)*u(2:n+1,j-1)-u(2:n+1,j-2);
    plot(x,u(:,j));
    axis([0,1,-1,1]);
    MM=getframe;
end

%plot 
mesh(u);