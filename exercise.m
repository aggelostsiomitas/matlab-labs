clear;
clc;
A=zeros(100,100); % μαυρο τετραγωνακι 
imshow(A);
A=ones(100,100); %ασπρο τετραγωνακι 
imshow(A);
A=rand(100,100);%
imshow(A);

N=99;
K=zeros(N,N);
v=ones(1,N);
r=(N+1)/2;
c=(N+1)/2;
xx=ones(N,N);
yy=ones(N,N);

for i=1:N
    xx(i,:)=i*xx(i,:);
    yy(:,i)=i*yy(:,i);
end

imshow( (xx-c).^2 +(yy-c).^2<=r^2);
%το παρακατω φταιχνει ενα πινακα με rand  και με υο getframe φαινεται να
%αλλαξει το χρωμα απο ασπρο σε μαυρο επανηλημενα 
 %for i=1:1000
%A=rand(N,N);
%imshow(A.*((xx-c).^2 +(yy-c).^2<=r^2));
%M=getframe;
 %end

