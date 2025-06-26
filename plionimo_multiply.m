clear;
clc;
a=input("put the first polionimo");
b=input("put the second polionimo");

n1=length(a);
n2=length(b);

disp(0:n1-1);
disp(a);

disp(0:n2-1);
 disp(b);

%γινομενο πινακα α με b με την χρηση της συναρτησης poliadd
if n1>n2 
    c=b(1)*a;
    for i=2:n2
        c=poliadd(c,b(i)*[zeros(1,i-1) a]);
    end
else
    c=a(1)*b;
  for i=2:n1
     c=poliadd(c,a(i)*[zeros(1,i-1) b]);
  end
end
disp([(1:length(c))' c']);
%ή θα εκανα b'*a

t=input("give a polonimo to find paragogo");
t1=length(t);

dt=t(2:t1);
dt=dt.*(1:t1-1);
disp(dt);

da=a;
for i=2:n1
    da=da(2:n1);
    da=da.*(1:n1-i+1);
    disp([(0:length(da)-1)' da']);
end






    