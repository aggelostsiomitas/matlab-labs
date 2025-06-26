clear;
clc;
x=0:0.001:20;
x1=0:0.3:20;
t=1:6;
f=@(x)exp(cos(pi*x))-6*log(1+abs(sin(x)));
y=f(x);

local_max=islocalmax(y);
local_min=islocalmin(y);

hold on;
plot(x,f(x)),title('Διάγραμμα της f(x)');
xlabel('x');
ylabel('f(x)');
plot(x(local_max),y(local_max),'or');
plot(x(local_min), y(local_min),'or');
plot(x(local_max),y(local_max),'k--');
plot(x(local_min), y(local_min),'k--');
hold off;

[X1,T]=meshgrid(x1,t);
F = exp(cos(pi * X1 .* T)) - 6 * T .* log(1 + abs(sin(X1)));
figure;
surf(F);



