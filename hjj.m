function [H] = hjj(X,j)
%function that returns the enthaply value for a stage(j) of a given matrix
%X; returns a single number/value
%used in the calculation of the derivative matrices(A,B,C)
r=9.5; B=62;
xv=X(1:5,:); xl=X(7:11,:);
sv=sum(xv); sl=sum(xl);
if j>=2&&j<=18
H=-sl(j-1)*hl(xl(:,j-1)./sl(j-1),X(6,j-1))-sv(j+1)*hv_dept(xv(:,j+1)./sv(j+1),X(6,j+1))+sl(j)*hl(xl(:,j)./sl(j),X(6,j))+sv(j)*hv_dept(xv(:,j)./sv(j),X(6,j));        
elseif j==1
H= sl(:,1)-(r/(r+1))*sv(:,1);
elseif j==19
H=sl(:,19)-B;
end

