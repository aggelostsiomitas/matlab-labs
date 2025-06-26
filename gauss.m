function x = gauss(A,b)
AE=[A b];
n=length(A);

for i=1:n-1
    for j=n:-1:i+1
        p=-AE(j,i)/A(i,i);
        for k=1:n+1
            A(j,k)=A(j,k)+p*A(i,k); % αυτο το αποτελεσμα ειναι ο πινακας L
        end
    end
    x=AE(1:n,1:n)\AE(:,n+1);% αυτο το αποτελεσμα ειναι ο πινακας U
end
%συνολικα L*U=A
