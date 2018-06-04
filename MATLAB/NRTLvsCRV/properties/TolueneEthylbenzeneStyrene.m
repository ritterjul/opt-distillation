BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
C = [2;3;4];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);
    
%alpha = [3.083969465648855 1.343511450381679 1]; %paper -> already very good
%alpha = [3.066252447253957 1.369532253533125 1.020242664518052]; %fit to curves
%alpha = [3.068081324066286 1.367558080445769 1.016693439678843]; %fit to points -> no change

alpha = alpha/alpha(3);