clc;
clear;
Vi=10;
a=pi/4;
Tmax=3;
nt=100;
dt=Tmax/nt;
g=9.8;
Vx=zeros(1,nt+1);
Vy=zeros(1,nt+1);
x=zeros(1,nt+1);
y=zeros(1,nt+1);
x(1)=0;
y(1)=0;
Vx(1)=Vi*cos(a);
Vy(1)=Vi*sin(a);

i=2;
while(y(i-1)>=0)
    Vx(i)=Vx(i-1);
    Vy(i)=Vy(i-1)-g*dt;
    x(i)=x(i-1) + Vx(i-1)*dt;
    y(i)=y(i-1) + Vy(i-1)*dt-1/2*g*(dt)^2;
    i=i+1;
end
plot(x(1:i-1),y(1:i-1));