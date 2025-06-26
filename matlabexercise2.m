clc;
clear;
Vi=10;
a=pi/4;
Tmax=1;
nt=100;
dt=Tmax/nt;
g=9.8;
xi=0;
yi=0;
Vxi=Vi*cos(a);
Vyi=Vi*sin(a);
Vx=zeros(1,nt+1);
Vy=zeros(1,nt+1);
x=zeros(1,nt+1);
y=zeros(1,nt+1);
%άμα ειχα απο το 1 στο nt το γραφημα δεn θα ξεκινουσε ακριβως απο το 0
%γιατι ο κωδικας ξεκινα οχι για t=0 αλλα για t=dt για αυτο βαζουμε n+1
for i=1:nt
    Vx(i)=Vxi;
    Vy(i)=Vyi-g*dt;
    x(i)=xi + Vxi*dt;
    y(i)=yi + Vyi*dt-1/2*g*(dt)^2;
    Vxi=Vx(i);
    Vyi=Vy(i);
    xi=x(i);
    yi=y(i);
end
disp(xi);
disp(yi);
plot(x,y);

    
 
