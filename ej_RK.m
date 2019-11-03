function [E] = ej_RK(X,j)
%function that returns equilibrium equations value for a X matrix
c=5; n=19; eff=0.85;
Tc=[513 508 510 562 537]; %Kelvin
Pc=[81 48 47.50 48.9 53.2868]; %bar
R=8.314;
P=101325; %in N/m2
bi=0.08664*R.*Tc./Pc;
E=zeros(c,1);
%antoine equation parameters:
a=[5.20409 4.42448 4.20364 4.72583 4.20772];
b=[1581.341 1312.253 1164.426 1660.652 1233.129];
C=[-33.50 -32.445 -52.69 -1.461 -40.953];
xv=X(1:5,:); xl=X(7:11,:); 
P=1.013;
sv=sum(xv); sl=sum(xl);
%liquid molar volumes for all components:
vl1=[40.73; 74.05; 79.84;89.41;80.67]; %in cm^3/mol
vl=vl1.*10^-3; %in m^3/kmol
%wilson parameters
param=[0.0 -2.1501*10^2 -9.8420*10^1 2.1613*10^2 -3.7225*10^2;
6.6429*10^2 0.0 1.6126*10^2 -1.6790*10^2 -3.1310*10^2;
8.34583*10^2 -6.5210*10^1 0.0 2.0054*10^2 -4.3141*10^2;
1.67946*10^3 4.949199*10^2 1.01700*10^1 0.0 -2.0850*10^2;
1.70257*10^3 -7.3150*10^1 7.91000*10^1 1.4844*10^2 0.0];
param=transpose(param);
R=1.987; %in cal/mol*Kelvin
A=zeros(5,5);
for k=1:5
    for i=1:5
        volrat(i,k)=vl(k)/vl(i);
    end
end
if (j>=2 && j<=18)
    A=volrat.*exp(-param./(R.*X(6,j)));
    x=(xl(:,j)./sl(j));
    sum_xA=A*x;
    ai=0.42748*R^2.*(Tc.^2.5)./(Pc.*X(6,j)^0.5);
    amix=0;
    for q=1:5
        for w=1:5
        amix=amix+(xv(q,j)/sv(j))*(xv(w,j)/sv(j))*(ai(q)*ai(w))^0.5;
        end
    end
    bmix=0;
    for q=1:5
        bmix=bmix+xv(q,j)*bi(q)/sv(j);
    end
    for i=1:5
        Psat=10^(a(i)-b(i)/(X(6,j)+C(i)));
        big_sum=0;
        for k=1:5
            big_sum=big_sum+(xl(k,j)./sl(j))*A(k,i)/sum_xA(k);
        end
        yi=-log(sum_xA(i)) + 1 - big_sum;
        p1= P*X(6,j)^0.5;
        p2= -1*R*X(6,j)^1.5;
        p3= amix-bmix*R*X(6,j)^1.5-(bmix^2)*(P)*(X(6,j)^0.5);
        p4= -1*amix*bmix;
        p= [p1 p2 p3 p4]; %p1,p1,p3,p4 are the coefficients of the cubic equation
        r= roots(p);
        v=max(r); %taking maximum of the roots as the gas phase molar volume
        phi = 2.718^((bi(i)*(P*v/(R*X(6,j))-1))/bmix-log((v-bmix)*P/(R*X(6,j)))+log(1+bmix/v)*(amix*bi(i)/bmix-2*(sum(xv(i,j)/sv(j))))/(bmix*R*X(6,j)^1.5));
        E(i,1)=-1*xv(i,j)+(Psat/P)*exp(yi)*xl(i,j)*sv(j)/(phi*sl(j))+(1-eff)*sv(j)*xv(i,j+1)/sv(j+1);
    end
end
if j==19 || j==1
    A=volrat.*exp(-param./(R.*X(6,j)));
    x=(xl(:,j)./sl(j));
    sum_xA=A*x;
    ai=0.42748*R^2.*(Tc.^2.5)./(Pc.*X(6,j)^0.5);
    amix=0;
    for q=1:5
        for w=1:5
        amix=amix+(xv(q,j)/sv(j))*(xv(w,j)/sv(j))*(ai(q)*ai(w))^0.5;
        end
    end
    bmix=0;
    for q=1:5
        bmix=bmix+xv(q,j)*bi(q)/sv(j);
    end
    for i=1:5
        Psat=10^(a(i)-b(i)/(X(6,j)+C(i)));
        big_sum=0;
        for k=1:5
            big_sum=big_sum+(xl(k,j)./sl(j))*A(k,i)/sum_xA(k);
        end
        yi=-log(sum_xA(i)) + 1 - big_sum;
        p1= P*X(6,j)^0.5;
        p2= -1*R*X(6,j)^1.5;
        p3= amix-bmix*R*X(6,j)^1.5-(bmix^2)*(P)*(X(6,j)^0.5);
        p4= -1*amix*bmix;
        p= [p1 p2 p3 p4]; %p1,p1,p3,p4 are the coefficients of the cubic equation
        r= roots(p);
        v=max(r); %taking maximum of the roots as the gas phase molar volume
        phi = 2.718^((bi(i)*(P*v/(R*X(6,j))-1))/bmix-log((v-bmix)*P/(R*X(6,j)))+log(1+bmix/v)*(amix*bi(i)/bmix-2*(sum(xv(i,j)/sv(j))))/(bmix*R*X(6,j)^1.5));
        E(i,1)=-1*xv(i,j)+(Psat/P)*exp(yi)*xl(i,j)*sv(j)/(phi*sl(j));
end
end

