function [P] = langrange(x,xi,yi)
n=length(xi);
n2=length(x);
P=zeros(1,n2);
    for j = 1:n
        Li = ones(1,n2);
        for k= [1:j-1 j+1:n]
                Li = Li.* (x - xi(k)) / (xi(j) - xi(k));
        end
        P = P + Li * yi(j);
    end
end