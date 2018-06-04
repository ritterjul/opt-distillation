clear;
%choose mixture
BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
letter = 'A';
%EthanolIso_1PropanolIso_1Butanol;
%letter = 'B';
%HexaneHeptaneOctaneNonaneDecane;
%letter = 'C';

%choose subset of components(optional)
C = [3,4,5];
number = strcat(num2str(C(1)),num2str(C(2)),num2str(C(3)));
cPr = cPr(:,C);
cG = cG(:,C,C);
cH = cH(:,C);
ch = ch(:,C);
alpha = alpha(C);
alpha = alpha/alpha(3);

%%
y_rig_all = zeros(3,36);
x_rig_all = zeros(3,36);
for j=1:8
    for k=1:(9-j)
        f = [j*0.1 k*0.1 1-(j+k)*0.1];
        y_rig_all(:,sum(9-(1:(j-1)))+k) = y_rig(f,Pr,cG,cPr);
        x_rig_all(:,sum(9-(1:(j-1)))+k) = x_rig(f,Pr,cG,cPr);
    end
end
diff = sum(sum(sum((y_rig_all-y_CRV_all(alpha)).^2+(x_rig_all-x_CRV_all(alpha)).^2)));

alpha_fit = fminunc(@(alpha)sum(sum(sum((y_rig_all-y_CRV_all(alpha)).^2+(x_rig_all-x_CRV_all(alpha)).^2))),alpha);
diff_fit = sum(sum(sum((y_rig_all-y_CRV_all(alpha_fit)).^2+(x_rig_all-x_CRV_all(alpha_fit)).^2)));
alpha_fit_reg = alpha_fit/alpha_fit(3);

save(strcat(letter,'_',number),'alpha','diff','alpha_fit','diff_fit','alpha_fit_reg');

function y = y_CRV_all(alpha)
y = zeros(3,36);
for j=1:8
    for k=1:(9-j)
        f = [j*0.1 k*0.1 1-(j+k)*0.1];
        y(:,sum(9-(1:(j-1)))+k)=y_CRV(f,alpha);
    end
end
end

function x = x_CRV_all(alpha)
x = zeros(3,36);
for j=1:8
    for k=1:(9-j)
        f = [j*0.1 k*0.1 1-(j+k)*0.1];
        x(:,sum(9-(1:(j-1)))+k)=x_CRV(f,alpha);
    end
end
end
