function [derh] = derh(X,j,s)
%function that returns the derivative of Hj equations with respect to vij,Tj
%and lij
c=5;
derh=zeros(1,11);
Tij=X(6,:);
for i=1:c
    %taking numerical derivative of Hj wrt to vij using central difference
    %method
    delv=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+delv;
    X22=X;
    X22(i,s)=X(i,s)-delv;
    derh(1,i)=(hjj(X11,j)-hjj(X22,j))/(2*delv);  %central difference
    
end
for i=6
    %taking numerical derivative of Hj wrt to Tj using central difference
    %method
    delt=Tij(1,s)/10^4;
    X11=X;
    X11(6,s)=X(6,s)+delt;
    X22=X;
    X22(6,s)=X(6,s)-delt;
    derh(1,6)=(hjj(X11,j)-hjj(X22,j))/(2*delt);  %central difference
end
for i=7:11
    %taking numerical derivative of Hj wrt to lij using central difference
    %method
    dell=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+dell;
    X22=X;
    X22(i,s)=X(i,s)-dell;
    derh(1,i)=(hjj(X11,j)-hjj(X22,j))/(2*dell);  %central difference
end

