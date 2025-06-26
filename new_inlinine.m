%function [c]=new_inlinine(s1,s2)
%c=inline([s1 '+' s2]);
function [c]=new_inlinine(s1,s2) % μπορουσα να βαλω στην παρενθεση το d
%function [c]=new_inlinine(s1,s2,d)
% και το τρεχω ετσι:new_inlinine('x.^2','e^x','+')
d=input("put a character");
c=inline([s1 d s2]);
end

