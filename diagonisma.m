%%
clc;
clear;

t=-2:0.01:9;
r1=(t+1).*heaviside(t+1);
r2=(t-2).*heaviside(t-2);
r3=(t-3).*heaviside(t-3);
r4=(t-5).*heaviside(t-5);
r5=(t-6).*heaviside(t-6);

r=r1-2*r2+r3-2*r4+2*r5;
plot(t,r);
xlim([-2.5 9.5]);
ylim([-0.5 3.5]);
%%
clc;
clear;
a=[32 0 -2];
b=[1 -2];
n=0:19;
x=[1 2 -1 zeros(1,17)];
y=filter(b,a,x);
stem(n,y);

x2=gauspuls(n);
y2=filter(b,a,x2);
figure;
stem(n,y2);

%%
clc;
clear;

num3=[4 0];
den3=[2 2 4];
H3=tf(num3,den3);

num4=[10 2 0];
den4=[1 2 5];
H4=tf(num4,den4);

H34=parallel(H3,H4);

num1=10;
den1=[2 1 4];
H1=tf(num1,den1);

num2=[4 0 1];
den2=[1 2 1];
H2=tf(num2,den2);

H12=series(H1,H2);

H=series(H12,H34)

poloi=pole(H)
zeros=zero(H)
pzmap(H);
xlim([-2 2]);

%οι πολοι βρισκονται αριστερλα του φανταστικού άξονα  αρα ειναι ευσταθες 

t=0:0.01:20;
h=impulse(H,t);
figure;
plot(t,h);

s=step(H,t);
figure;
plot(t,s);


x3=2*t.*exp(-10.*t-2).*heaviside(t-4);
y3=lsim(H,x3,t);
figure;
plot(t,y3);