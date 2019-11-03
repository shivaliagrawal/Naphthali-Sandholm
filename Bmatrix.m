function [B] = Bmatrix(X,f,j)
%function that returns the B matrix 
%in the dF/dX tridiagonal matrix
H=derh(X,j,j);
M=derm(X,f,j,j);
E=dere(X,j,j);
B=[H; M; E];
end
