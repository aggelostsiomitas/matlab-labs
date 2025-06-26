clear;
clc;
xi=0:0.1:1;
n=length(xi);
yi=[0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.75 0.7 0.71 0.64];
x=linspace(0,1,100);
n2=length(x);
P=langrange(x,xi,yi);

plot(xi,yi,'or');hold
plot(x,P); hold off

h=x(2)-(x(1));
I=sum(P(1:n2-1))*h;
disp(I);

   