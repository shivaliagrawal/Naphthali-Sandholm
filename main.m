%Shivali Agrawal 160651
%Assignment No.6

c=5; %no. of components
n=19; %no. of stages inside the column
X1=zeros(2*c+1,n); %initial guess matrix
reflux=9.5;

%1-Methanol
%2-Acetone
%3-Methyl Acetate
%4-Benzene
%5-Chloroform

%Antoine parameters for all the components:
P=1.013; %in bar
A=[5.20409 4.42448 4.20364 4.72583 4.20772];
B=[1581.341 1312.253 1164.426 1660.652 1233.129];
C=[-33.50 -32.445 -52.69 -1.461 -40.953];

Tlow=329.4107; %boiling point of maximum volatile component
Thigh=353.2776; %boiling point of minimum volatile component

for i=1:n
    X1(c+1,i)=(Tlow*(n-i)+Thigh*(i-1))/(n-1); %interpolating the temperature of all the stages between Tmax at the lowest stage and Tmin at the topmost stage
end
Tflash=(Tlow*(n-7)+Thigh*(7-1))/(n-1); %taking temperature of the feed plate as the flash temperature
Psat=zeros(1,c);
k=zeros(1,c);
%using Rachford Rice for the initial guess of V and L and assuming Raoults law for calculating ki
for i=1:c
    Psat(i)=10^(A(i)-B(i)/(Tflash+C(i)));
    k(i)=Psat(i)/P;
end
phiguess=0;
phinew=0.9;
z=[0.15; 0.40; 0.05; 0.20; 0.20]; %feed mole fractions
while (abs(phinew-phiguess)>0.001)
    phiguess=phinew;
    f=z(1)*(k(1)-1)/(1+phiguess*(k(1)-1))+z(2)*(k(2)-1)/(1+phiguess*(k(2)-1))+z(3)*(k(3)-1)/(1+phiguess*(k(3)-1))+z(4)*(k(4)-1)/(1+phiguess*(k(4)-1))+z(5)*(k(5)-1)/(1+phiguess*(k(5)-1));
    fder=-1*z(1)*(k(1)-1)^2/(1+phiguess*(k(1)-1))^2-1*z(2)*(k(2)-1)^2/(1+phiguess*(k(2)-1))^2-1*z(3)*(k(3)-1)^2/(1+phiguess*(k(3)-1))^2-1*z(4)*(k(4)-1)^2/(1+phiguess*(k(4)-1))^2-1*z(5)*(k(5)-1)^2/(1+phiguess*(k(5)-1))^2;
    phinew=phiguess-f/fder;
end

feed=100;
V=phinew*feed;
L=feed-V;

x=zeros(1,c);
y=zeros(1,c);

for i=1:c
    x(i)=z(i)/((k(i)-1)*phinew+1);
    y(i)=k(i)*x(i);
end
l=x.*L;
v=y.*V;

for i=1:5
    X1(i,:)=v(i);
    X1(i+6,:)=l(i);
end

X=horzcat(X1(:,1)',X1(:,2)');
for i=3:19
    X= horzcat(X, X1(:,i)');
end

F=zeros(2*c+1,n);
f=zeros(c,n);
for i=1:5
    f(i,7)=z(i)*feed;
end
s=zeros(1,n);
s(1)=reflux;

lambda=hv_dept(z,330);
tau=2; epsilon=1;
iter=0;
while tau>epsilon
    disp('iteration')
   disp(iter);
   disp('tau')
    disp(tau)
    disp('epsilon')
     disp(epsilon)
H=hj(X1); %enthalpy balance equations
M=Material(X1,f); %material balance equations
E=Eij_RK(X1); %equilibrium equations

F2=horzcat(H',M');
F=horzcat(F2,E');
F=F';

F1=horzcat(F(:,1)',F(:,2)');
for i=3:19
    F1= horzcat(F1, F(:,i)');
end

N=19;
size=2*c+1;
dfdx=zeros(size*N,size*N);

%creating the df/dx matrix as follows:
for j=1:n-1
    dfdx(1+size*(j-1):size*j,j*size+1:size*(j+1))=Cmatrix(X1,f,j);
end
for j=1:n
    dfdx(1+size*(j-1):size*j,1+size*(j-1):size*j)=Bmatrix(X1,f,j);
end
for j=2:n
    dfdx(1+size*(j-1):size*j,1+size*(j-2):size*(j-1))=Amatrix(X1,f,j);
end

%updating the X matrix using Newton raphson algorithm
Xnew=X'-inv(dfdx)*F1';
X=Xnew';
for j=1:N
    X1(:,j)=X((2*c+1)*(j-1)+1:(2*c+1)*(j));
end
 
epsilon=size*N*(F1*F1')*10^(-10); %updating epsilon
%calculating tau:
 for j=1:N
     F1(size*(j-1)+1)= F1(size*(j-1)+1)/lambda;
 end
 tau=F1*F1';
 iter=iter+1;
end
%X1(6,:)=X1(6,:)+n+1;

for i=1:19
    Nstage(i)=i;
end
figure(1)
plot(Nstage,X1(6,:));
title('Block temperature profile')
ylabel('Temperature (Kelvin)')
xlabel('Stage')
legend('Temperature (Kelvin)')
xv=X1(1:5,:); xl=X1(7:11,:); 
sv=sum(xv); sl=sum(xl);
ys=zeros(5,19);
for j=1:19
for i=1:5
ys(i,j)=xv(i,j)/sv(j);
end
end
xs=zeros(5,19);
for j=1:19
for i=1:5
xs(i,j)=xl(i,j)/sl(j);
end
end

figure(2)
plot(Nstage,xs)
title('Liquid composition profiles')
ylabel('Mole fraction')
xlabel('Stage')
legend('liquid mole fraction of methanol','liquid mole fraction of acetone','liquid mole fraction of methyl acetate','liquid mole fraction of benzene','liquid mole fraction of chloroform')
figure(3)
plot(Nstage,ys)
title('Vapour composition profiles')
ylabel('Mole fraction')
xlabel('Stage')
legend('vapour mole fraction of methanol','vapour mole fraction of acetone','vapour mole fraction of methyl acetate','vapour mole fraction of benzene','vapour mole fraction of chloroform')














