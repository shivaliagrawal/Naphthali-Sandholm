function [dere] = dere(X,j,s) 
%function that returns the derivative of Eij_RK equations with respect to vij,Tj
%and lij
c=5;
dere=zeros(c,2*c+1);
Tij=X(6,:);
for k=1:c
for i=1:5 %taking numerical derivative of Eij_RK wrt to vij using central difference method
    delv=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+delv;
    X22=X;
    X22(i,s)=X(i,s)-delv;
    Mj1=ej_RK(X11,j);
    Mj2=ej_RK(X22,j);
    dere(k,i)=(Mj1(k,1)-Mj2(k,1))/(2*delv); %central difference
end
for i=6
    %taking numerical derivative of Eij_RK wrt to Tj using central difference
    %method
    delt=Tij(1,s)/10^4;
    X11=X;
    X11(6,s)=X(6,s)+delt;
    X22=X;
    X22(6,s)=X(6,s)-delt;
    Mj1=ej_RK(X11,j);
    Mj2=ej_RK(X22,j);
    dere(k,6)=(Mj1(k,1)-Mj2(k,1))/(2*delt); %central difference
end
for i=7:11
    %taking numerical derivative of Eij_RK wrt to lij using central difference
    %method
    dell=X(i,s)/10^4;
    X11=X;
    X11(i,s)=X(i,s)+dell;
    X22=X;
    X22(i,s)=X(i,s)-dell;
    Mj1=ej_RK(X11,j);
    Mj2=ej_RK(X22,j);
    dere(k,i)=(Mj1(k,1)-Mj2(k,1))/(2*dell); %central difference
end
end
end
