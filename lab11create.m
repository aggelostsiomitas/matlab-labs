clear;
clc;
n=10000;
for i=1:n
student(i).AEM=round(55000+(59000-55000)*rand);
student(i).physics=1+round(rand*9);%τυχαιοι αριθμοι απο το 1 μεχρι απο το 9
%a,b για τυχαιους αριθμους απο α στο β α+(β-α)*rand
end
c=0;
for i=1:n
    if student(i).physics>=5
        c=c+1;
    end
end
disp(c);