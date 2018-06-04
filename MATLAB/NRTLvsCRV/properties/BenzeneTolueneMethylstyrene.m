BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
C = [1;2;5];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);
    
%alpha = [10.5 4.04 1.0]; %paper -> bad
%alpha = [10.455042319755405 4.239214217195130 0.663140443807239]; %fit to curves
%alpha = [10.355591026136850 4.493047191313384 0.657758742859536]; %fit to points -> still bad

alpha = alpha/alpha(3);