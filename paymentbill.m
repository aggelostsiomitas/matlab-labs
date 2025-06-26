clear;
clc;
C=[126 501 689 1562 1781];
 n=length(C);
 cost=10*ones(1,n);
 
for i=1:n
  if C(i)<=500
    cost(i) =10+0.02*C(i);  
  elseif C(i)<=1000
    cost(i)=10+ 10+0.05*(C(i)-500);    
  elseif C(i)>1000
    cost(i)=10+35+0.1*(C(i)-1000);    
  end
end

%disp(['cost:' num2str(cost)]); %το κανω string
disp('cost:');
disp(cost');





    
