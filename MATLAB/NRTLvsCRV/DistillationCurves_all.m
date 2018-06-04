
clear;
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
            %if exist('alpha','var')
                %alpha = alpha(C);
                %alpha = alpha/alpha(3);
            %end
            load(strcat('results\fit\',letter,'_',number,'.mat'),'alpha_fit_reg');
            alpha = alpha_fit_reg;
            %figure
            %ternplot(1,0,0)
            %hold on
            diff = 0;
            for j=1:3
                for k = 1:(4-j)
                    rig = zeros(11,3);
                    rig(6,:) = [j*0.2 k*0.2 1-(j+k)*0.2];
                    for i=5:-1:1
                        rig(i,:) = x_rig(rig(i+1,:),Pr,cG,cPr);
                    end
                    for i=7:11
                        rig(i,:) = y_rig(rig(i-1,:),Pr,cG,cPr);
                    end
                    %ternplot(rig(:,1),rig(:,2),rig(:,3),'r-*')
                    %if exist('alpha','var')
                        CRV = zeros(11,3);
                        CRV(6,:) = [j*0.2 k*0.2 1-(j+k)*0.2];
                        for i=5:-1:1
                            CRV(i,:) = x_CRV(CRV(i+1,:),alpha);
                        end
                        for i=7:11
                            CRV(i,:) = y_CRV(CRV(i-1,:),alpha);
                        end
                        ternplot(CRV(:,1),CRV(:,2),CRV(:,3),'b-*')
                        diff = diff + sum(sum((rig-CRV).^2));
                   %end
                end
            end
            %hold off
            %savefig(strcat(letter,'_',number,'.fig'));
            save(strcat(letter,'_',number),'diff');
        end
    end
end
close all;
