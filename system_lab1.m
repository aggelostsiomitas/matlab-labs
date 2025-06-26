%%
clear;
clc;
t=-3:0.01:3;
u1=heaviside(t);
u2=heaviside(t-2);
x=t.^3.*cos(10*pi*t).*(u1-u2);
plot(t,x);

%%
x=-10:4:10;
y=x.^2;
plot(x,y);
x2=-10:1:10;
y2=x2.^2;
hold on;
plot(x2,y2);

%%
x=-10:4:10;
y=x.^2;
plot(x,y);
x2=-10:1:10;
y2=x2.^2;
hold on;
stem(x2,y2); % κανει ενα διαγραμμα το οποιο ενωνει με γραμουλες τα σημεια μου 

%%
t1=linspace(0,1,1000);
t2=linspace(1,6,1000);
y1=t1;
y2=1./t2;
t=[t1 t2];
y=[y1 y2];
plot(t,y);

%%
t=linspace(0,10,1000);
x=t.*exp(-t).*cos(2*pi*4*t);
plot(t,x);
xlabel('Time(s)');
ylabel('Amptitude');

%%
t=-7:0.01:7;
u=heaviside(t);
plot(t,u);
hold on
u2=heaviside(t-2);
plot(t,u2);
hold on
u3=heaviside(t+2);
plot(t,u3);

%%
t=-5:.1:10;
s=gauspuls(t); % value=1 for t=0
plot(t,s);

%%
t=-5:.1:10;
r=t.*heaviside(t);
plot(t,r);
%y=x equation so angle =45 degrees


%%
clc;
clear;
syms t 
u=heaviside(t);
int(u,t);

%%
clc;
clear;
t=-3:0.01:3;
u1=heaviside(t+2);
u2=heaviside(t-2);
p=u1-u2;
subplot(3,1,1);
plot(t,u1);
subplot(3,1,2);
plot(t,u2);
subplot(3,1,3);
plot(t,p);

%%
clc;
clear;
n=-3:3;
u=heaviside(n);
stem(n,u);

%%
n=-3:3;
s=gauspuls(n);
stem(n,s);

%%
n=-3:3;
r=n.*heaviside(n);
stem(n,r);

%%
n=-10:10;
x=(0.9.^n) .*exp(j*n);
subplot(4,1,1);
stem(n,real(x));
subplot(4,1,2);
stem(n,abs(x));
subplot(4,1,3);
stem(n,angle(x));
subplot(4,1,4);
stem(n,imag(x));

%%
t=-3:0.01:3;
u1=heaviside(t);
u2=heaviside(t-2);
x=t.^3.*cos(10*pi*t).*(u1-u2);
plot(t,x);

%%
clc;
clear;
t=-1:0.01:3;
r1=t.*heaviside(t);
r2=(t-1).*heaviside(t-1);
r3=(t-2).*heaviside(t-2);
r=r1-r2-r3;
plot(t,r);

%%
clc;
clear;

t=-1:0.01:10;
r1=t.*heaviside(t);
r2=(t-1).*heaviside(t-1);
r3=(t-2).*heaviside(t-2);
r4=(t-3).*heaviside(t-3);
r5=(t-5).*heaviside(t-5);
r=r1-r2+r3-2*r4+r5;
plot(t,r);
xlim([-2 10]);% puts limits in my plot this allows me to see it more clearly 
%or limit what i can see
ylim([-0.3 2.2]);

%%
clc;
clear;

t=-3:0.01:3;
x=-heaviside(t-1)+heaviside(t-2);

r1= (t+1).*heaviside(t+1);
r2=heaviside(t-2);
h=2*r1.*heaviside(2-t);

y=conv(x,h)*0.01;
t_plot=-6:0.01:6;
plot(t_plot,y);
ylim([-6 3]);
xlim([-8 8]);
%%
clc;
clear;

t = -3:0.01:3;
x = -heaviside(t-1) + heaviside(t-2);

% r1 = (t+1) .* heaviside(t+1);
% r2 = heaviside(t-2);
% h = 2 * r1 .* heaviside(2-t);
h = 2 * (t+1) .* heaviside(t+1) - 2 * (t+1) .* heaviside(t-2);

y = conv(x, h) * 0.01; % Υπολογισμός συνέλιξης

t_conv = linspace(min(t) + min(t), max(t) + max(t), length(y)); % Νέος άξονας χρόνου

plot(t_conv, y); % Σχεδίαση με σωστό μήκος
ylim([-6 3]);
xlabel('Χρόνος');
ylabel('Πλάτος');
title('Συνέλιξη x(t) * h(t)');
grid on;

%%
clc;
clear;
t = -2:0.01:2; % Ορισμός του χρονικού άξονα
x = (t >= -2 & t <= 0) .* 1 + (t > 0 & t <= 2) .* exp(-3*t); % Ορισμός του x(t)
h = heaviside(t) - heaviside(t - 3); % Κρουστική απόκριση h(t)

y=conv(x,h)*0.01;
t_plot=-4:0.01:4;
plot(t_plot,y);
xlim([-4.5 4.5]);
ylim([-2 3]);