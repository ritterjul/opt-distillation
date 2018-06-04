HexaneHeptaneOctaneNonaneDecane;
C = [1;2;3];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);
    
%alpha =[2.78 1.78 1]; %paper
%alpha = [3.269290711498254 1.495035250776431 0.734117694981484]; %fit curves -> better but not good
%alpha = [3.141416777369321 1.433193889509823 0.702974366345384]; %fit points -> better but not good

alpha = alpha/alpha(3);


