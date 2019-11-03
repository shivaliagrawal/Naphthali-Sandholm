function [derm] = derm(X,f,j,s) 
%function that returns the derivative of Mij equations with respect to vij,Tj
%and lij
c=5;
derm=zeros(c,2*c+1);
Tij=X(6,:);
for k=1:5
for i=1:5
    %taking numerical derivative of Mij wrt to vij using central difference
    %method
    delv=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+delv;
    X22=X;
    X22(i,s)=X(i,s)-delv;
    Mj1=Mj(X11,f,j);
    Mj2=Mj(X22,f,j);
    derm(k,i)=(Mj1(k,1)-Mj2(k,1))/(2*delv); %central difference
end
for i=6
    %taking numerical derivative of Mij wrt to Tj using central difference
    %method
    delt=Tij(1,s)/10^4;
    X11=X;
    X11(6,s)=X(6,s)+delt;
    X22=X;
    X22(6,s)=X(6,s)-delt;
    Mj1=Mj(X11,f,j);
    Mj2=Mj(X22,f,j);
    derm(k,6)=(Mj1(k,1)-Mj2(k,1))/(2*delt); %central difference
end
for i=7:11
    %taking numerical derivative of Mij wrt to lij using central difference
    %method
    dell=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+dell;
    X22=X;
    X22(i,s)=X(i,s)-dell;
    Mj1=Mj(X11,f,s);
    Mj2=Mj(X22,f,s);
    derm(k,i)=(Mj1(k,1)-Mj2(k,1))/(2*dell); %central difference
end
end
end

