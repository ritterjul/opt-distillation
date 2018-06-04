function [x] = x_CRV(y,alpha)
x=(1./alpha).*y/((1./alpha)*y');
end
