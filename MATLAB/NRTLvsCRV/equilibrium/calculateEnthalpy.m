function [H,h] = calculateEnthalpy(T,Tref,cH,ch)
Href = (cH(1,:)*Tref + cH(2,:).*cH(3,:)./tanh(cH(3,:)/Tref) + cH(4,:).*cH(5,:).*tanh(cH(5,:)/Tref))/3.6e6;
H = (cH(1,:)*T + cH(2,:).*cH(3,:)./tanh(cH(3,:)/T) + cH(4,:).*cH(5,:).*tanh(cH(5,:)/T))/3.6e6-Href;
h =  H - (ch(1,:).*(1-T./ch(2,:)).^(ch(3,:) + T./ch(2,:).*(ch(4,:) + T./ch(2,:).* ch(5,:))))/3.6e6;
end

