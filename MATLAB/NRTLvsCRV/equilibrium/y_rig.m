function [y] = y_rig(x,Pr,cG,cPr)
T = fzero(@(T)(calculatePressure(T,cPr).*calculateGamma(x,T,cG)*x'-Pr),350);
y = calculatePressure(T,cPr)/Pr.*calculateGamma(x,T,cG).*x;
end

