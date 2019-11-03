function [ M ] = Mj(X,f,j)
%function that returns material balance (C*1 matrix) for a particular j value
%used in making the derivative matrices (A,B,C)
c=5;
reflux=9.5;
M=zeros(c,1);
if j>=2&&j<=18
    for i=1:5
    M(i,1)=-1*X(c+1+i,j-1)-1*X(i,j+1)-f(i,j)+X(c+1+i,j)+X(i,j);
    end
elseif j==1
for i=1:5
    M(i,1)=-1*X(i,2)+X(c+1+i,1)*(1+1/reflux)+X(i,1);
end
elseif j==19
    for i=1:5
    M(i,1)=-1*X(c+1+i,18)+X(c+1+i,19)+X(i,19);
    end
end
end


