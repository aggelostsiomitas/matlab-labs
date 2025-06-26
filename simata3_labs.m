%%
clear;
clc;

syms t s;
u=heaviside(t);
h=t*exp(-t)*u;
H=laplace(h,s);

%%
clear;
clc;

num=[1 0 0];
den=[1 3 1 ];
Hs=tf(num,den);

%%
clear;
clc;
syms s t;
H=s/(s^2+4);
h=ilaplace(H,t);
ezplot(h,[0 20]);

%%
clear;
clc;

den=[1 3 6 4];
poles=roots(den);
plot(poles, '+')

num=[3 0 5];
midenika=roots(num);
hold on
plot(midenika,'o')
xlim([-2 2]);

%%
%δευτερος τροπος αντι για το παραπανω 
clear;
clc;

den=[1 3 6 4];
num=[3 0 5];
H=tf(num,den);
poloi=pole(H);
midenika=zero(H);
pzmap(H);
xlim([-2 2]);

%%
clear;
clc;

num1 = 1;
den1 = [500 0 0];
num2 = [1 1];
den2 = [1 2];
H1 = tf(num1, den1);
H2 = tf(num2, den2);
H = series(H1, H2);
H_par=parallel(H1,H2);
H_FD=feedback(H1,H2);

%%
clear;
clc;

num=10;
den=[1 2 10];
t=0:0.01:10;
x=cos(2*pi*t);
y=lsim(num,den,x,t);
plot(t,y)

%%
clear;
clc;

num=10;
den=[1 2 10];
t=0:0.01:10;
x1=heaviside(t);
y1=lsim(num,den,x1,t);
plot(t,y1)
figure
x2=gauspuls(t);
y2=lsim(num,den,x2,t);
plot(t,y2)

%%
%αλλος τροπος αντι για το παραπανω 
clear;
clc;

num=10;
den=[1 2 10];
t=0:0.01:10;

y1=step(num,den,t);
plot(t,y1)
figure

y2=impulse(num,den,t);
plot(t,y2)

%%
clear;
clc;

syms n z
h=2^n;
H=ztrans(h,z);

%%
clear;
clc;

num=[2 1];
den=[1 3 2];
midenika=roots(num);
poloi=roots(den);
zplane(midenika,poloi)

%%
clear;
clc;

num=[3 -1.4 0.15];
den=[1 -0.7 0.15 -0.025];
midenika=roots(num);
poloi=roots(den);
zplane(midenika,poloi)

%second solution
H=tf(num,den,0.1)
pzmap(H)

%third solution 
metra=abs(poloi)

%%
clear;
clc;

num=[.1 .1];
den=[1 -1.5 .7];
[y1 n1]=stepz(num,den)
stem(n1,y1)

[y2 n2]=impz(num,den)
figure
stem(n2,y2)

%%
clear;
clc;

num = [.1 -.1];
den = [1 -1.5 0.7];
n = 0:50;
x = gauspuls(n);
y = dlsim(num, den, x);
stem(n, y)

%%
clear;
clc;
num=[.1 .1];
den=[1 -1.5 0.7];
n=0:50;
x=(-1).^n;
y=dlsim(num,den,x);
stem(n,y)

%%
clear;
clc;

num1 = [1 +5];
num2 = [1 -2];
num = conv(num1, num2);
den1 = [4 j];
den2 = [4 -j];
den = conv(den1, den2);
H = tf(num, den);

%%
clc;
clear;

num1=[2 0 1];
den1=[1 3 3 1];
H1=tf(num1,den1);

n1=[1 1];
n2=[1 2];
num2=conv(n1,n2);

d1=[1 2j];
d2=[1 -2j];
d3=[1 3];
den2=conv(d1,d2);
den2_final=conv(den2,d3);
%θα μπορουσα να γραψω den2=conv(conv([1 2j],[1-2j]),[1 3]);
H2=tf(num2,den2_final);

H_total=tf(H1,H2);  
poloi=pole(H_total);
pzmap(H_total)

%%
clc;
clear;

num1=[1 0];
den1=[1 4];
H1=tf(num1,den1);

num2=[1 0];
den2=[1 0 4];
H2=tf(num2,den2);

num3=[1 0 0];
den3=[1 0 4];
H3=tf(num3,den3);

H12=parallel(H1,H2);
H=series(H12,H3);

%%
clc;
clear;

num1=conv([1 2],[1 -3]);
den1=[1 9 26 24];
H=tf(num1,den1);
t=0:0.1:3;
x=gauspuls(t);
y=lsim(num1,den1,x,t);
plot(t,y)


%%
clc;
clear;
syms s
H=(s+2)*(s-3)/(s^3+9*s^2+26*s+24);
h=ilaplace(H);
ezplot(h,[0 3])
ylim([-0.4 1])

%%
clc;
clear;

t=0:60;
x=exp(-0.1*t).*cos(t);
num=[1 0 -1];
den=[1 2 3 4 ];
Y=lsim(num,den,x,t);
plot(t,Y)

%%
clc;
clear;

num=[8 10 -6];
den=[1 2 -1 -2];
midenika=roots(num);
poloi=roots(den);
zplane(midenika,poloi)
metra=abs(poloi)

%%
clc;
clear;

num=[1 0 0];
den=[1 -5/6 1/6];
n=0:13;
x=gauspuls(n);
y=dlsim(num,den,x);
stem(n,y)

%%
clc;
clear;

num=0.0798*ones(1,4);
den=[1 -1.556 1.272 -0.398];
midenika=roots(num);
poloi=roots(den);
zplane(midenika,poloi)
metra=abs(poloi)

