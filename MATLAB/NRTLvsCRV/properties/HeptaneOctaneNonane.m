HexaneHeptaneOctaneNonaneDecane;
C = [2;3;4];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);
    
%alpha =[2.535211267605634 1.478873239436620 1]; %paper -> okay
%alpha = [2.708399906766324 1.362735689073346 0.730566195842380]; %fit curves -> very good
%alpha = [2.708399906766324 1.362735689073346 0.730566195842380]; %fit points -> very good

alpha = alpha/alpha(3);


