function [C] = Cmatrix(X,f,j)
%function that returns the C 
%matrix in the dF/dX tridiagonal matrix
H=derh(X,j,j+1);
M=derm(X,f,j,j+1);
E=dere(X,j,j+1);
C1=[H; M; E];
W=zeros(11,5);
C=horzcat(C1(:,1:6),W);
end
