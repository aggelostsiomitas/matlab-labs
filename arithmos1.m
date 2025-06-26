clear;
clc;
n=6;
counter=1;
while(n>1)
    if mod(n,2)==0
        n=n/2;
    else
        n=3*n+1;
    end
    counter=counter+1;
end
disp(counter);