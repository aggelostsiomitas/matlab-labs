%%
clc;
clear;
step=0.01;
t=0:step:2;
x=heaviside(t);
t1=0:step:1;
h=1-t1;
y=conv(x,h)*step;
t_plot=0:step:3;
plot(t_plot,y);
%%
clear;
clc;
syms t
h = exp(t.^2);
int(abs(h), t, -inf, inf)


%%
clc;
clear;
t=-10:0.01:10;
x=heaviside(t+10)-heaviside(t-10);
h=exp(-t.^2);
y=conv(x,h)*0.01;
t_plot=-20:0.01:20;
plot(t_plot,y);

%%
clc;
clear;
syms t
h=t.^2;
int(abs(h), t, -inf, inf)

%%
clc;
clear;
t=-10:0.01:10;
h=t.^2;

x=heaviside(t+10)-heaviside(t-10);
y=conv(x,h)*0.01;
t_plot=-20:0.01:20;
plot(t_plot,y);


%%
clear;
clc;
n=0:3; 
h=[2 4 3 1];
x=gauspuls(n);
y=conv(x,h);

n_stem=0:6;
stem(n_stem,y);


%%
clc;
clear;
syms n
h = 1/(2^n);
symsum(abs(h), n, 0, inf)

%%
clc;
clear;

N=10;
a=[1 -0.5 .1];
b=[0.1 -0.5 1];
x=[1 2 -1 zeros(1,7)]; %προσθετω 7 μεδενικα γιατι θελω 10 σημεια αλλα εχω μονο 3
y=filter(b,a,x);
n=0:N-1;
stem(n,y);

%%
clc;
clear;

a=[1 -0.5 .1];
b=[0.1 -0.5 1];
n=0:9;  
x=gauspuls(n);
y=filter(b,a,x);
stem(n,y);

%%
clc;
clear;
a=[1 -1];
b=[.2 .1];
n=0:100;
x=heaviside(n);
y=filter(b,a,x);
stem(n,y);

title('Step response of the system')

%%
clc;
clear;

t=0:0.01:10;
h=t;
x=0.8.^t;

y=conv(x,h)*0.01;
t_plot=0:0.01:20;
plot(t_plot,y);

%%
clc;
clear;
t=0:0.01:5;
x=heaviside(t-2)-heaviside(t-4);
h=heaviside(t)-heaviside(t-2);

y=conv(x,h)*0.01;
t_plot=0:0.01:10;
plot(t_plot,y);



%%
clc;
clear;
n=0:10;
h=0.7.^n;
x=heaviside(n)-heaviside(n-4);

y=conv(x,h);
n_stem=0:20;
stem(n_stem,y);


%%
clc;
clear;
syms t
h1=exp(-3*t);
h2=t.*exp(-2*t);
h=h1+h2;
int(abs(h),t,0,inf)


%%
clc;
clear;
%αν ηθελα συνολικη αποκριση σοτ προηγουμενο 
t=0:0.01:10;
x=heaviside(t);
h=exp(-3*t) + t.*exp(-2*t);
y=conv(x,h)*0.01;
t_plot=0:0.01:20;
plot(t_plot,y);

%%
clc;
clear;
n=0:20;
a=[1 0 -0.8];
b=[1 -0.5];
x=gauspuls(n);
y=filter(b,a,x);

x2=heaviside(n);
y2=filter(b,a,x2);
stem(n,y);
figure;
stem(n,y2)

%%
