clc;
clear;
U=10;
a=pi/4;
t=0:1:5;
g=9.8;
Ux=U*cos(a);
Uy=U*sin(a);
x=Ux*t;
y=Uy*t-1/2*g*t.^2;
disp(x);
disp(y);
plot(x,y);

