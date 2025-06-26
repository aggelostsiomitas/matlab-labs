clear;
clc;

x = 0:0.001:20;
f=@(x)exp(cos(pi*x))-6*log(1+abs(sin(x)));
y = f(x);

figure;
plot(x, y, 'b'); 
hold on;

[local_max, local_min] = find_local(x, y);

plot(local_max, f(local_max), 'ro'); 
plot(local_min, f(local_min), 'ro');

title('Γράφημα της f(x) με τοπικά μέγιστα και ελάχιστα');
xlabel('x');
ylabel('f(x)');
grid on;

plot(local_max, f(local_max),'k--');
plot(local_min, f(local_min),'k--');
hold off;

t = 1:0.1:6; 
x1=0:0.5:20;
[X, T] = meshgrid(x1, t);

Z = exp(cos(pi * X .* T)) - 6 * T .* log(1 + abs(sin(X)));


figure;
mesh(X, T, Z); 
title('Γράφημα της συνάρτησης f(x,t)');
xlabel('x');
ylabel('t');
zlabel('f(x,t)');
grid on;
