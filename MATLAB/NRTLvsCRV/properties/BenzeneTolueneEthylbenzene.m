BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
C = [1;2;3];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);

alpha = alpha(C);
    
%alpha = [5.9659 2.2955 1]; %paper -> already very good
%alpha = [6.020186996568187 2.175807565621761 0.951761991203437]; %fit to curves
%alpha = [6.020436092277860 2.179078139177632 0.948444785426739]; %fit to points -> no real change

alpha = alpha/alpha(3);