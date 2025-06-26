clear;
clc;
theta_start=5*pi/360;
nmax=4;
nmin=1;
l=5;
h=(nmax-nmin)/(l-1);
n=nmax-(0:l-1)*h;
theta=theta_start;
x=0;
y=0;
figure;
hold on;
xlabel("x");
ylabel("y");
grid on;
for i=1:l-1
  theta_critical = asin(n(i+1) / n(i));
   if theta>theta_critical
       disp("ΕΧΟΥΜΕ ΟΛΙΚΗ ΑΝΑΚΛΑΣΗ")
       break;
   else
   theta=asin(n(i)/n(i+1)*sin(theta));
   newX=x+cos(theta);
   newY=y+sin(theta);
   plot(x,y,newX,newY,'ob--');
    plot([x, newX], [y, newY], 'r-');

   end
   x=newX;
   y=newY;
end
hold off;


    
    