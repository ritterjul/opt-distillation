EthanolIso_1PropanolIso_1Butanol;
C = [2;3;5];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);

%alpha = [3.6 2.1 1]; %paper -> already good
%alpha = [3.644585566235112 2.042722823355697 0.962954923104727]; %fit curves
%alpha = [3.640447258463868 2.050335168492005 0.958547030353904]; %fit points -> slight change, maybe better?

alpha = alpha/alpha(3);