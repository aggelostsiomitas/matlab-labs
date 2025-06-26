clear;
clc;
a=16;
num=[];

while (a~=0)
    num=[num2str(mod(a,2)) num];
     a=floor(a/2);
end
disp(num);







