function pS = calculatePressure(T,cPr)
pS = exp(cPr(1,:) + cPr(2,:)/T + cPr(3,:)*log(T) + cPr(4,:).*T.^cPr(5,:));    
end