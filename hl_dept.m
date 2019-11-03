function [h] = hl_dept( x, Tf )
%function that returns the liquid molar enthaplies at the given temperature
% comp 1-methanol, 2-acetone, 3-methyl acetate, 4-benzene, 5-chloroform
%heat capacities in the liquid phase
%in J/kmol*kelvin
cpl{1}=@(T) (105800-362.23.*T+0.9379.*T.^2);
cpl{2}=@(T) (135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4);
cpl{3}=@(T) (61260+270.9.*T);
cpl{4}=@(T) (162940-344.94.*T+0.85562.*T.^2);
cpl{5}=@(T) (124850-166.34.*T+0.43209.*T.^2);
h=0;
for i=1:5
    %T=330K is taken as the reference temperature
    h=h+x(i)*integral(cpl{i},330, Tf); %in J/kmol
    h=h/1000;
end
end