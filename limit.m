clear;
clc;
x=1;
N=20;
e=exp(1);
for i=1:N
   c=(1+x)^(1/x);
  disp([c abs(c-e)]);
  x=x/10;
end

format long






