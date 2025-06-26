clear;
clc;
N=100000;
c=zeros(1,N);
for i=1:N
    c(i)=functionarithmos(i);
end
figure;
plot(c,'k');

figure;
histogram(c,20);

[v,i]=max(c);%θα επιστρεψει την τιμη και την θεση του μεγιστου του c
disp([v,i]);