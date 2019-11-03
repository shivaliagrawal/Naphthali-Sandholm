function [A] = Amatrix(X,f,j)
%function that returns the A 
%matrix in the dF/dX tridiagonal matrix
c=5;
H=derh(X,j,j-1);
A1=zeros(c,c+1);
I=-eye(c);
A2=zeros(c,2*c+1);
A3=horzcat(A1,I);
A=vertcat(H,A3,A2);
end
