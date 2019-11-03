function [enthalpy] = VapourEnthalpy( y,Temp )

%Enthalpy of vapour calculation ideal case
 c = 5; % No. of components

Latent_Heat = [35.21 29.1 30.32 30.72 29.24]*10^6; % J/Kmol
Tsat = [337.8 329 330.1 353.1 334.2]; % in Kelvin

% Liquid Heat Capacity 

cpl{1}=@(T) (105800-362.23.*T+0.9379.*T.^2);
cpl{2}=@(T) (135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4);
cpl{3}=@(T) (61260+270.9.*T);
cpl{4}=@(T) (162940-344.94.*T+0.85562.*T.^2);
cpl{5}=@(T) (124850-166.34.*T+0.43209.*T.^2);

% Vapor Heat Capacity 

cpv{1}=@(T) (0.39252*10^5+0.879.*10^5.*((1.9165.*(10^3)./T)./sinh(1.9165.*(10^3)./T)).^2+0.53654.*10^5.*((896.7./T)./sinh(896.7./T)).^2);
cpv{2}=@(T) (0.5704*10^5+1.632.*10^5.*((1.607.*(10^3)./T)./sinh(1.607.*10^3./T)).^2+0.968.*10^5.*((731.5./T)./sinh(731.5./T)).^2);
cpv{3}=@(T) (0.555*10^5+1.782.*10^5.*((1.26.*(10^3)./T)./sinh(1.26.*10^3./T)).^2+0.853.*10^5.*((562./T)./sinh(562./T)).^2);
cpv{4}=@(T) (0.44767*10^5+2.3085.*10^5.*((1.4792.*(10^3)./T)./sinh(1.4792.*10^3./T)).^2+1.6836.*10^5.*((677.66./T)./sinh(677.66./T)).^2);
cpv{5}=@(T) (0.3942*10^5+0.6573.*10^5.*((0.928.*(10^3)./T)./sinh(0.928.*10^3./T)).^2+0.493.*10^5.*((399.6./T)./sinh(399.6./T)).^2);


enthalpy = 0;

for i = 1:5
    enthalpy= enthalpy + y(i) * (integral(cpl{i},330, Tsat(i)) + Latent_Heat(i) + integral(cpv{i},Tsat(i), Temp)/1000); %J/Kmol
    
end


% adding non-ideal part using Departure Functions

Pc = [81 48 47.50 48.9 53.28]*0.986923; %in atm critical Pressure
Tc = [513 508 510 562 537];             %in Kelvin critical temperature

R=0.0821;             %in  ltr-atm/mol-k

% pure component parameter apure , bpure
apure=0.42748*(R^2.*(Tc).^2.5./Pc);
bpure=0.08664*R.*Tc./Pc;

A = zeros(c); % mixture parameter

for i = 1:c
    for j=1:c
       A(i,j)=sqrt(apure(i)*apure(j));
    end
end

a=0;

for i = 1:c
    for j=1:c
       a = a+y(i)*y(j)*A(i,j);
    end
end

b = sum(y.*bpure');

% ### Equations are taken from Seader and Henley Separation Process Chapter2

P=1;    %in atm
A=a*P/(R^2*Temp^2);
B=b*P/(R*Temp);

% for equation we need to calculate v as we  calculated in ASSIGNMENT2

eq=[P*Temp^0.5 -1*R*Temp^1.5 a-b*R*Temp^1.5-(b^2)*(P)*(Temp^0.5) -1*a*b];
root=roots(eq);
v=max(root);
Zv=(P*v)/(R*Temp);
Departure = (R*Temp*(Zv-1-3*A*log(1+B/Zv)/(2*B))); 
enthalpy=enthalpy+Departure;

end

