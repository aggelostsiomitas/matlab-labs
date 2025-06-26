%%
clear;
clc;

%κουτακια
m=4;

%number of potitions 
n=(m-1)^2;

%domain size,diastima
L=[-1,1];

%mesh step 
h=(L(2)-L(1))/m;%((b-a)/m)

%function
G=1;
sigma=0.1;
f=@(x,y) 4*pi*G*(1/sqrt(pi*sigma))*exp(-(x.^2+y.^2)/(2*sigma^2));

%coordinates
x=linspace(L(1)+h,L(2)-h,m-1);
y=linspace(L(1)+h,L(2)-h,m-1);

[xx,yy]=ndgrid(x,y);
%αντι για το παραπανω θα μπορουσ ανα κανω αυτο:
xx=reshape(xx,n,1);
yy=reshape(yy,n,1);

%coetrix matrix 
A=zeros(n,n);
for i=1:n
    A(i,i)=-4;
end

%diagonals before and after the main diagonal
for i=1:n-1
    A(i,i+1)=1;
    A(i+1,i)=1;
   % if mod(i,m-1)==0
    %A(i,i+1)=0;
    %A(i+1,i)=0;
   % end
end
for i=1:m-2
    A(i*(m-1),i*(m-1)+1)=0;
    A(i*(m-1)+1,i*(m-1))=0;
end

%further diagonals 
for i=1:n-(m-1)
    A(i,m+i-1)=1;
    A(m+i-1,i)=1;
end

b=h^2*f(xx,yy);
x=A\b;

%plot the solution 
mesh(reshape(xx,m-1,m-1),reshape(yy,m-1,m-1),reshape(x,m-1,m-1));
%%