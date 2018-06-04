clear;
%choose mixture
%BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
%EthanolIso_1PropanolIso_1Butanol;
%AcetoneChloroformBenzene;
HexaneHeptaneOctaneNonaneDecane;
%MethanolEthanolWater;
%AcetoneEthylacetateEthanol;

%choose subset of components(optional)
C = [3;4;5];
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
if exist('alpha','var')
    alpha = alpha(C);
    alpha = alpha/alpha(3);
end
%%
figure
ternplot(1,0,0)
hold on
for j=1:8
    for k=1:(9-j)
        f = [j*0.1 k*0.1 1-(j+k)*0.1];
        rig = y_rig(f,Pr,cG,cPr);
        ternplot(rig(1),rig(2),rig(3),'r*')
        if exist('alpha','var')
            CRV = y_CRV(f,alpha);
            ternplot(CRV(1),CRV(2),CRV(3),'b*')          
        end
        rig = x_rig(f,Pr,cG,cPr);
        ternplot(rig(1),rig(2),rig(3),'r*')
        if exist('alpha','var')
            CRV = x_CRV(f,alpha);
            ternplot(CRV(1),CRV(2),CRV(3),'b*')           
        end
    end
end
hold off

