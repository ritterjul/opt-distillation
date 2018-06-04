function [x] = x_rig(y,Pr,cG,cPr)
options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 1000;
z = fsolve(@(z)([calculatePressure(z(4),cPr).*calculateGamma([z(1) z(2) z(3)],z(4),cG).*[z(1) z(2) z(3)]-Pr*y, sum(z(1:3))-1]),[y 350],options);
x = z(1:3);
end


 