% x- CONFIDENTIAL

Pr = x;

cPr = [ x x x x x;
		x x x x x;
        x x x x x;
        x x x x x;
        x x x x x];
			
cG = zeros(3,5,5);
cG(2,:,:) = x;
cG(:,1,4) = [x; x; x];
cG(:,2,5) = [x; x; x];
cG(:,4,1) = [x; x; x];
cG(:,5,2) = [x; x; x];
  
Tref = 298;

cH = [ x x x x x;
		x x x x x;
        x x x x x;
        x x x x x;
        x x x x x];
    
ch = [ x x x x x;
		x x x x x;
        x x x x x;
        x x x x x;
        x x x x x];
    
alpha = [10.5 4.04 1.76 1.31 1.0];