function gamma = calculateGamma(x,T,cG)
tau = squeeze(cG(2,:,:)+cG(3,:,:)/T);
P = exp(-squeeze(cG(1,:,:)).*tau);
S1 = (P.*tau)'*x';
S2 = P'*x';
gamma = exp(S1./S2 + P.*(tau - (S1./S2)')*(x'./S2))';
end
