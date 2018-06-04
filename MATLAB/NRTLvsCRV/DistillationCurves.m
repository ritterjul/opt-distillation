clear;
%choose mixture
%BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
%EthanolIso_1PropanolIso_1Butanol;
HexaneHeptaneOctaneNonaneDecane;
%AcetoneChloroformBenzene;
%MethanolEthanolWater;
%AcetoneEthylacetateEthanol;

%choose subset of components(optional)
C = [2;3;4];
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
diff = 0;
for j=1:3
    for k=1:(4-j)
        rig = zeros(11,3);
        rig(6,:) = [j*0.2 k*0.2 1-(j+k)*0.2];
        for i=5:-1:1
            rig(i,:) = x_rig(rig(i+1,:),Pr,cG,cPr);
        end
        for i=7:11
            rig(i,:) = y_rig(rig(i-1,:),Pr,cG,cPr);
        end
        ternplot(rig(:,1),rig(:,2),rig(:,3),'r-*')
        if exist('alpha','var')
            CRV = zeros(11,3);
            CRV(6,:) = [j*0.2 k*0.2 1-(j+k)*0.2];
            for i=5:-1:1
                CRV(i,:) = x_CRV(CRV(i+1,:),alpha);
            end
            for i=7:11
                CRV(i,:) = y_CRV(CRV(i-1,:),alpha);
            end
            ternplot(CRV(:,1),CRV(:,2),CRV(:,3),'b-*')
        end
        diff = diff + sum(sum((CRV-rig).^2));
    end
end
hold off