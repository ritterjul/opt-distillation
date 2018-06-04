EthanolIso_1PropanolIso_1Butanol;
C = [1;3;5];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);

%alpha = [4.1 2.1 1.0]; %paper -> good
%alpha = [4.157484519056275 2.013801375372557 0.948015483869656]; %fit curves
%alpha = [4.156018168250838 2.016675072600953 0.945158775320002]; %fit points -> slight change, maybe better?

alpha = alpha/alpha(3);