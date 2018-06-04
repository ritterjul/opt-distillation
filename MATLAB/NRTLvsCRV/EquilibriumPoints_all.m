
for a = 1:3
    for b = a+1:4
        for c = b+1:5
            %choose mixture
            %BenzeneTolueneEthylbenzeneStyreneMethylstyrene;
            %letter = 'A';
            %EthanolIso_1PropanolIso_1Butanol;
            %letter = 'B';
            HexaneHeptaneOctaneNonaneDecane;
            letter = 'C';
            C = [a;b;c];
            number = strcat(num2str(a),num2str(b),num2str(c));
            cPr = cPr(:,C);
            cG = cG(:,C,C);
            cH = cH(:,C);
            ch = ch(:,C);
            if exist('alpha','var')
                alpha = alpha(C);
                alpha = alpha/alpha(3);
            end
            figure
            ternplot(1,0,0)
            hold on
            diff = 0;
            for j=1:8
                for k=1:(9-j)
                    f = [j*0.1 k*0.1 1-(j+k)*0.1];
                    rig = y_rig(f,Pr,cG,cPr);
                    ternplot(rig(1),rig(2),rig(3),'r*')
                    if exist('alpha','var')
                        CRV = y_CRV(f,alpha);
                        ternplot(CRV(1),CRV(2),CRV(3),'b*')
                        diff = diff + sum((rig-CRV).^2);
                    end
                    rig = x_rig(f,Pr,cG,cPr);
                    ternplot(rig(1),rig(2),rig(3),'r*')
                    if exist('alpha','var')
                        CRV = x_CRV(f,alpha);
                        ternplot(CRV(1),CRV(2),CRV(3),'b*')
                        diff = diff + sum((rig-CRV).^2);
                    end
                end
            end
            hold off
            savefig(strcat(letter,'_',number,'.fig'));
            save(strcat(letter,'_',number),'diff');
        end
    end
end
close all;
