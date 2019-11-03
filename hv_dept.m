function [h] = hv_dept( y, Tf )
%function that returns vapour enthalpies 
% comp 1-methanol, 2-acetone, 3-methyl acetate, 4-benzene, 5-chloroform

Tc=[513 508 510 562 537]; %Kelvin
Pc=[81 48 47.50 48.9 53.2868]; %bar
R=8.314;
P=101325; %in N/m2
bi=0.08664*R.*Tc./Pc;
ai=0.42748*R^2.*(Tc.^2.5)./(Pc.*Tf^0.5);
a=0;
for i=1:5
    for j=1:5
        a=a+y(i)*y(j)*(ai(i)*ai(j))^0.5;
    end
end
b=0;
for i=1:5
    b=b+y(i)*bi(i);
end
A=a*P/(R^2*Tf^2);
B=b*P/(R*Tf);
%heat capacities in the liquid phase
cpl{1}=@(T) (105800-362.23.*T+0.9379.*T.^2);
cpl{2}=@(T) (135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4);
cpl{3}=@(T) (61260+270.9.*T);
cpl{4}=@(T) (162940-344.94.*T+0.85562.*T.^2);
cpl{5}=@(T) (124850-166.34.*T+0.43209.*T.^2);

%heat capacities in the vapour phase
cp{1}=@(T) (0.39252*10^5+0.879.*10^5.*((1.9165.*(10^3)./T)./sinh(1.9165.*(10^3)./T)).^2+0.53654.*10^5.*((896.7./T)./cosh(896.7./T)).^2);
cp{2}=@(T) (0.5704*10^5+1.632.*10^5.*((1.607.*(10^3)./T)./sinh(1.607.*10^3./T)).^2+0.968.*10^5.*((731.5./T)./cosh(731.5./T)).^2);
cp{3}=@(T) (0.555*10^5+1.782.*10^5.*((1.26.*(10^3)./T)./sinh(1.26.*10^3./T)).^2+0.853.*10^5.*((562./T)./cosh(562./T)).^2);
cp{4}=@(T) (0.44767*10^5+2.3085.*10^5.*((1.4792.*(10^3)./T)./sinh(1.4792.*10^3./T)).^2+1.6836.*10^5.*((677.66./T)./cosh(677.66./T)).^2);
cp{5}=@(T) (0.3942*10^5+0.6573.*10^5.*((0.928.*(10^3)./T)./sinh(0.928.*10^3./T)).^2+0.493.*10^5.*((399.6./T)./cosh(399.6./T)).^2);

Hvap=[71.8 31.3 35.2 30.72 29.4]; %kJ/mol
Hvap=10^6.*Hvap; % J/Kmol
Tsat=[64.7 56 57.1 80.1 61.2]; % degree C
Tsat=Tsat+273.15; % Kelvin
h=0;
  %enthalpies calculated are in J/Kmol
    %T=330K is taken as the reference temperature
    %enthaply in vapour phase is calculated as the enthalpy of liquid in
    %going from T=330 to Tboil plus the enthalpy of vapourisation at Tboil
    %plus the enthalpy of vapour in going from Tboil to given temperature
for i=1:5
    h=h+y(i)*(integral(cpl{i},330, Tsat(i)) + Hvap(i) + integral(cp{i},Tsat(i), Tf));
    h=h/1000;
end
        p1= P*Tf^0.5;
        p2= -1*R*Tf^1.5;
        p3= a-b*R*Tf^1.5-(b^2)*(P)*(Tf^0.5);
        p4= -1*a*b;
        p= [p1 p2 p3 p4]; %p1,p1,p3,p4 are the coefficients of the cubic equation
        r= roots(p);
        v=max(r); %taking maximum of the roots as the gas phase molar volume
        Zv=P*v/(R*Tf);
        h=h+R*Tf*(Zv-1-3*A*log(1+B/Zv)/(2*B));
end

