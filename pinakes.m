clear;
clc;
n=3;
m=5;
k=4;
A=rand(n,m);
B=rand(m,k);
N=rand(m,1);
disp(A);
disp(B);
disp(N);
tic;
for i=1:n
    for j=1:k
     % s=0;
     % for p=1:m
      % s=s+A(i,p).*B(p,j);
      % end
      %  C(i,j)=s;
      % τον παραπανω κωδικα μπορω να τον απλοποιησω βγαζοντας την εξτρα for
      %μπορω να γραψω αυτο 
       %C(i,j)=sum(A(i,1:m).*Β(1:m,j)');
       %μπορωομως να γραψω κατευθειαν 
      % C(i, j) = A(i,:).* B(:,j); το * υποδηλωνει εσωτερικο γινομενο
      % αρα η τελεια μπορει να παραληφθει¨
      C(i, j) = A(i,:) * B(:,j);
    end
end
t=toc;
disp(['Elapsed Time:' num2str(t)]);%% γινομενο πινακων
%ακομα πιο γρηγορο εδω χρησιμοποιω πινακα επι δυανυσμα
for p=1:k
    t=zeros(n,1);
   for j=1:n
    t=t+A(:,j)*B(j,p);
    end
end

%τελικη γρηγορη λυση 
for p=1:m
    t=zeros(n,1);
    for j=1:m
        t=t+A(:,j)*B(j,p);
    end
    C(:,p)=t;
end
disp(['GFLOPS:' num2str(n*m*k/(10^9*t))]);


