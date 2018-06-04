#rigorous model of superstructure for two columns 

set C; #set of components

set I = {1,2}; #set of columns

param NS {I} integer, >= 1; #numbers of stages
set S {i in I} = 1..NS[i]; #sets of stages;
set Sm {i in I} = 2..NS[i]-1; #sets of middle stages
param fn integer; #feed stage, if fixed
param ccn integer; #stage for connecting stream from condenser
param crn integer; #stage for connecting stream from reboiler
param psn integer; #side product stage

param Fall >= 0; #feed flow
param xF {C} >= 0, <= 1; #feed concentrations
var F {Sm[1]} >= 0; 
subject to Feed: sum{n in Sm[1]} F[n] = Fall;
check: sum{c in C} xF[c] = 1;

var V {i in I, S[i]} >= 0; #vapor mass flow between stages
var y {i in I, S[i], C} >= 0, <= 1; #vapor concentrations between stages
subject to ClosureY {i in I, n in S[i]}: 1 = sum{c in C} y[i,n,c];

var L {i in I, S[i]} >= 0; #liquid mass flow between stages
var x {i in I, S[i], C} >= 0, <= 1; #liquid concentrations between stages
subject to ClosureX {i in I, n in S[i]}: 1 = sum{c in C} x[i,n,c];

var CC {Sm[2]} >= 0, <= Fall; #connecting mass flow from before condenser of first column
#same concentration as y[1,NS[1],c], same temperature, but liquid
var CR {Sm[2]} >= 0, <= Fall; #connecting mass flow from before reboiler of first column
#same concentration as x[1,1,c], same temperature and phase
var PC {I} >= 0, <= Fall; #product mass flow after condenser
#same concentration as y[i,NS[i],c], same temperature, but liquid
var PR {I} >= 0, <= Fall; #product mass flow before reboiler
#concentration: x[i,1,c], same temperature and phase
var PS {Sm[2]} >= 0, <= Fall; #product mass flow at side stream of second column
#same concentration as x[2,n,c], same temperature and phase 

var T {i in I, S[i]} >= 275, <= 500; #temperatures in stages
var TF >= 275, <= 500; #temperature feed

param Pr >= 0; #pressure in columns

var QC {I} >= 0, <= 10000*Fall; #heat duty condenser
var QR {I} >= 0, <= 10000*Fall; #heat duty reboiler

param cH {1..5, C}; #model constants
param Tref; #reference temperature
param Href {c in C} = (cH[1,c]*Tref + cH[2,c]*cH[3,c]/tanh(cH[3,c]/Tref) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/Tref)); #reference enthalphy
var H {i in I, n in S[i], c in C} = (cH[1,c]*T[i,n] + cH[2,c]*cH[3,c]/tanh(cH[3,c]/T[i,n]) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/T[i,n])-Href[c])/3.6e6; #enthalphy of vapor
var HF {c in C} = (cH[1,c]*TF + cH[2,c]*cH[3,c]/tanh(cH[3,c]/TF) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/TF)-Href[c])/3.6e6; #enthalphy of vapor

param ch {1..5, C}; #model constants
var h {i in I, n in S[i], c in C} =  H[i,n,c] - (ch[1,c]*(1-T[i,n]/ch[2,c])^(ch[3,c] + T[i,n]/ch[2,c]*(ch[4,c] + T[i,n]/ch[2,c]* ch[5,c])))/3.6e6; #enthalphy of liquid
var hF {c in C} = HF[c] - (ch[1,c]*(1-TF/ch[2,c])^(ch[3,c] + TF/ch[2,c]*(ch[4,c] + TF/ch[2,c]* ch[5,c])))/3.6e6; #enthalphy of liquid

param cGamma {1..3, C, C}; #model constants
var tauGamma {i in I, n in S[i], c in C, d in C} = (cGamma[2,c,d] + cGamma[3,c,d]/T[i,n]);
var PGamma {i in I, n in S[i], c in C, d in C} = exp(-cGamma[1,c,d]*tauGamma[i,n,c,d]);
var SGamma1 {i in I, n in S[i], c in C} = sum{d in C}(x[i,n,d]*PGamma[i,n,d,c]*tauGamma[i,n,d,c]);
var SGamma2 {i in I, n in S[i], c in C} = sum{d in C}(x[i,n,d]*PGamma[i,n,d,c]);
var gamma {i in I, n in S[i], c in C} = exp(SGamma1[i,n,c]/SGamma2[i,n,c] + sum{d in C}(x[i,n,d]*PGamma[i,n,c,d]*(tauGamma[i,n,c,d]-SGamma1[i,n,d]/SGamma2[i,n,d])/SGamma2[i,n,d])); #activity coefficient 

param cPr {1..6, C}; #model constants
var PrS {i in I, n in S[i], c in C} = exp(cPr[1,c] + cPr[2,c]/T[i,n] + cPr[3,c]*log(T[i,n]) + cPr[4,c]*T[i,n]^cPr[5,c]); #vapor pressure
var PrSF {c in C} = exp(cPr[1,c] + cPr[2,c]/TF + cPr[3,c]*log(TF) + cPr[4,c]*TF^cPr[5,c]);
subject to BoilingLiquidFeed: sum{c in C} xF[c]*PrSF[c] = Pr;

subject to ComponentBalanceTopStage1 {c in C}:
V[1,NS[1]-1]*y[1,NS[1]-1,c] = (PC[1] + sum{n in Sm[2]} CC[n])*y[1,NS[1],c] + L[1,NS[1]]*x[1,NS[1],c] ;
subject to ComponentBalance1 {n in Sm[1], c in C}:
V[1,n-1]*y[1,n-1,c] + L[1,n+1]*x[1,n+1,c] + F[n]*xF[c] = V[1,n]*y[1,n,c] + L[1,n]*x[1,n,c];
subject to ComponentBalanceBottomStage1 {c in C}:
L[1,2]*x[1,2,c] = V[1,1]*y[1,1,c] + (PR[1] + sum{n in Sm[2]} CR[n])*x[1,1,c];
 
subject to ComponentBalanceTopStage2 {c in C}:
V[2,NS[2]-1]*y[2,NS[2]-1,c] = PC[2]*y[2,NS[2],c] + L[2,NS[2]]*x[2,NS[2],c] ;
subject to ComponentBalance2 {n in Sm[2], c in C}:
V[2,n-1]*y[2,n-1,c] + L[2,n+1]*x[2,n+1,c] + CC[n]*y[1,NS[1],c] + CR[n]*x[1,1,c]= V[2,n]*y[2,n,c] + L[2,n]*x[2,n,c] + PS[n]*x[2,n,c];
subject to ComponentBalanceBottomStage2 {c in C}:
L[2,2]*x[2,2,c] = V[2,1]*y[2,1,c] + PR[2]*x[2,1,c] ;

subject to EnergyBalanceTopStage1:
sum{c in C} (V[1,NS[1]-1]*y[1,NS[1]-1,c]*H[1,NS[1]-1,c] ) = sum{c in C}((PC[1] + sum{n in Sm[2]} CC[n])*y[1,NS[1],c]*h[1,NS[1],c] + L[1,NS[1]]*x[1,NS[1],c]*h[1,NS[1],c]) + QC[1];
subject to EnergyBalance1 {n in Sm[1]}:
sum{c in C} (V[1,n-1]*y[1,n-1,c]*H[1,n-1,c] + L[1,n+1]*x[1,n+1,c]*h[1,n+1,c] + F[n]*xF[c]*hF[c])= sum{c in C} (V[1,n]*y[1,n,c]*H[1,n,c] + L[1,n]*x[1,n,c]*h[1,n,c] );
subject to EnergyBalanceBottomStage1:
sum{c in C} (L[1,2]*x[1,2,c]*h[1,2,c]) + QR[1]= sum{c in C} (V[1,1]*y[1,1,c]*H[1,1,c] + (PR[1] + sum{n in Sm[2]} CR[n])*x[1,1,c]*h[1,1,c]);

subject to EnergyBalanceTopStage2:
sum{c in C} (V[2,NS[2]-1]*y[2,NS[2]-1,c]*H[2,NS[2]-1,c]) = sum{c in C}(PC[2]*y[2,NS[2],c]*h[2,NS[2],c] + L[2,NS[2]]*x[2,NS[2],c]*h[2,NS[2],c]) + QC[2];
subject to EnergyBalance2 {n in Sm[2]}:
sum{c in C} (V[2,n-1]*y[2,n-1,c]*H[2,n-1,c] + L[2,n+1]*x[2,n+1,c]*h[2,n+1,c] + CC[n]*y[1,NS[1],c]*h[1,NS[1],c] + CR[n]*x[1,1,c]*h[1,1,c])= sum{c in C} (V[2,n]*y[2,n,c]*H[2,n,c] + L[2,n]*x[2,n,c]*h[2,n,c] + PS[n]*x[2,n,c]*h[2,n,c]);
subject to EnergyBalanceBottomStage2:
sum{c in C} (L[2,2]*x[2,2,c]*h[2,2,c]) + QR[2]= sum{c in C} (V[2,1]*y[2,1,c]*H[2,1,c] + PR[2]*x[2,1,c]*h[2,1,c]);

subject to Thermodynamics {i in I, n in S[i], c in C}:
0 = Pr*y[i,n,c] - PrS[i,n,c]*x[i,n,c]*gamma[i,n,c];

var Vap = V[1,1] + V[2,1];
minimize VapourFlow: Vap;





