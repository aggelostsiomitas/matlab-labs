function [s] = rightfill(s,n)
if length(s)< n
  s(length(s):n)=(' ');
end
if length(s)>=n
    s=s(1:n);
end
end

% για διαγραφη s(n=1:length(s))=[];
