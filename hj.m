function [ H ] = hj(X)
%function that returns the enthaply values for all the stages for a given X
%matrix; uses the prevoiusly defined functions hv_dept and hl
n=19;
H=zeros(1,n); r=9.5; 
B=62; %in kmol/hr
xv=X(1:5,:); xl=X(7:11,:);
sv=sum(xv); sl=sum(xl);

for j=2:18 %Q for all the stages from 2 to 18 is zero
H(j)=-sl(j-1)*hl(xl(:,j-1)./sl(j-1),X(6,j-1))-sv(j+1)*hv_dept(xv(:,j+1)./sv(j+1),X(6,j+1))+sl(j)*hl(xl(:,j)./sl(j),X(6,j))+sv(j)*hv_dept(xv(:,j)./sv(j),X(6,j));        
end
H(1)= sl(:,1)-(r/(r+1))*sv(:,1); %H1 and H19 equations are replaced by these equations because Q1 and Q19 are not given
H(19)= sl(:,19)-B;       
end

