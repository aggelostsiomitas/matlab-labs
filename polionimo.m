a=input("put the first polionimo");
b=input("put the second polionimo");
% το πολυωνυμο θα εκτυπωνει αναποδα τις δυναμεις πχ αν εχω ενα 5 βαθμου,
% το x^5 +2x^3+x-1 θα εκτυπωνει 0 1 2 3 4 5 
%θα εκτυπωνει -1 1 0 2 0 1
l1=length(a);
disp(0:l1-1);
disp(a);

l2=length(b);
disp(0:l2-1);
disp(b);
c=[b zeros(1,l1-l2)]+[a zeros(1,l2-l1)];
disp(c)