function [c] = poliadd(a,b)

l1=length(a);
l2=length(b);
c=[b zeros(1,l1-l2)]+[a zeros(1,l2-l1)];
end

