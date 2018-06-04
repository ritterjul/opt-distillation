% x- CONFIDENTIAL

Pr = x;

cPr = [ x x x x x;
		x x x x x;
        x x x x x;
        x x x x x;
        x x x x x];
			
cG = zeros(3,5,5);
cG(2,:,:) = x;
cG(:,1,2) = [x; x; x];
cG(:,2,1) = [x; x; x];
  
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
    
alpha =[8.9 5.7 3.2 1.55 1.0];



