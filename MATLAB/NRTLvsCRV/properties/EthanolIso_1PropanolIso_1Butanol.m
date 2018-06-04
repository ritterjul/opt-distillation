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
cG(:,2,3) = [x; x; x];
cG(:,2,4) = [x; x; x];
cG(:,3,2) = [x; x; x];
cG(:,3,4) = [x; x; x];
cG(:,4,1) = [x; x; x];
cG(:,4,2) = [x; x; x];
cG(:,4,3) = [x; x; x];
cG(:,4,5) = [x; x; x];
cG(:,5,4) = [x; x; x];
  
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

alpha = [4.1 3.6 2.1 1.42 1.0];