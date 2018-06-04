#rigorous model of column

set C; #set of components

param NS integer, >= 1; #numbers of stages
set S = 1..NS; #sets of stages;
set Sm = 2..NS-1; #sets of middle stages
param fn integer; #feed stage, if fixed

param Fall >= 0; #feed flow
param xF {C} >= 0, <= 1; #feed concentrations
var F {Sm} >= 0; 
subject to Feed: sum{n in Sm} F[n] = Fall;
check: sum{c in C} xF[c] = 1;

var V {S} >= 0; #vapor mass flow between stages
var y {S, C} >= 0, <= 1; #vapor concentrations between stages
subject to ClosureY {n in S}: 1 = sum{c in C} y[n,c];

var L {S} >= 0; #liquid mass flow between stages
var x {S, C} >= 0, <= 1; #liquid concentrations between stages
subject to ClosureX {n in S}: 1 = sum{c in C} x[n,c]; 

var PC >= 0, <= Fall; #product mass flow after condenser
#same composition as y[NS,i], same temperature, but liquid
var PR >= 0, <= Fall; #product mass flow before reboiler
#same composition as x[1,i], same temperature, liquid

var T {S} >= 275, <= 500; #temperatures in stages
var TF >= 275, <= 500; #temperature of feed

param Pr >= 0; #pressure in column

var QC >= 0, <= 1000*Fall; #heat duty condenser
var QR >= 0, <= 1000*Fall; #heat duty reboiler

param cH {1..5, C}; #model constants
param Tref; #reference temperature
param Href {c in C} = (cH[1,c]*Tref + cH[2,c]*cH[3,c]/tanh(cH[3,c]/Tref) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/Tref)); #reference enthalphy
var H {n in S, c in C} = (cH[1,c]*T[n] + cH[2,c]*cH[3,c]/tanh(cH[3,c]/T[n]) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/T[n])-Href[c])/3.6e6; #enthalphy of vapor
var HF {c in C} = (cH[1,c]*TF + cH[2,c]*cH[3,c]/tanh(cH[3,c]/TF) + cH[4,c]*cH[5,c]*tanh(cH[5,c]/TF)-Href[c])/3.6e6; #enthalphy of vapor

param ch {1..5, C}; #model constants
var h {n in S, c in C} =  H[n,c] - (ch[1,c]*(1-T[n]/ch[2,c])^(ch[3,c] + T[n]/ch[2,c]*(ch[4,c] + T[n]/ch[2,c]* ch[5,c])))/3.6e6; #enthalphy of liquid
var hF {c in C} = HF[c] - (ch[1,c]*(1-TF/ch[2,c])^(ch[3,c] + TF/ch[2,c]*(ch[4,c] + TF/ch[2,c]* ch[5,c])))/3.6e6; #enthalphy of liquid

param cGamma {1..3, C, C}; #model constants
var tauGamma {n in S, c in C, d in C} = (cGamma[2,c,d] + cGamma[3,c,d]/T[n]);
var PGamma {n in S, c in C, d in C} = exp(-cGamma[1,c,d]*tauGamma[n,c,d]);
var SGamma1 {n in S, c in C} = sum{d in C}(x[n,d]*PGamma[n,d,c]*tauGamma[n,d,c]);
var SGamma2 {n in S, c in C} = sum{d in C}(x[n,d]*PGamma[n,d,c]);
var gamma {n in S, c in C} = exp(SGamma1[n,c]/SGamma2[n,c] + sum{d in C}(x[n,d]*PGamma[n,c,d]*(tauGamma[n,c,d]-SGamma1[n,d]/SGamma2[n,d])/SGamma2[n,d])); #activity coefficient 

param cPr {1..6, C}; #model constants
var PrS {n in S, c in C} = exp(cPr[1,c] + cPr[2,c]/T[n] + cPr[3,c]*log(T[n]) + cPr[4,c]*T[n]^cPr[5,c]); #vapor pressure
var PrSF {c in C} = exp(cPr[1,c] + cPr[2,c]/TF + cPr[3,c]*log(TF) + cPr[4,c]*TF^cPr[5,c]);
subject to BoilingLiquidFeed: sum{c in C} xF[c]*PrSF[c] = Pr;

subject to ComponentBalanceTopStage {c in C}:
V[NS-1]*y[NS-1,c] = PC*y[NS,c] + L[NS]*x[NS,c] ;
subject to ComponentBalance {n in Sm, c in C}:
V[n-1]*y[n-1,c] + L[n+1]*x[n+1,c] + F[n]*xF[c] = V[n]*y[n,c] + L[n]*x[n,c];
subject to ComponentBalanceBottomStage {c in C}:
L[2]*x[2,c] = V[1]*y[1,c] + PR*x[1,c];
 
subject to EnergyBalanceTopStage1:
sum{c in C} (V[NS-1]*y[NS-1,c]*H[NS-1,c]) = sum{c in C} (PC*y[NS,c]*h[NS,c] + L[NS]*x[NS,c]*h[NS,c]) + QC;
subject to EnergyBalance1 {n in Sm}:
sum{c in C} (V[n-1]*y[n-1,c]*H[n-1,c] + L[n+1]*x[n+1,c]*h[n+1,c] + F[n]*xF[c]*hF[c])= sum{c in C} (V[n]*y[n,c]*H[n,c] + L[n]*x[n,c]*h[n,c]);
subject to EnergyBalanceBottomStage1:
sum{c in C} (L[2]*x[2,c]*h[2,c]) + QR = sum{c in C} (V[1]*y[1,c]*H[1,c] + PR*x[1,c]*h[1,c]);


subject to Thermodynamics {n in S, c in C}:
0 = Pr*y[n,c] - PrS[n,c]*x[n,c]*gamma[n,c];

#minimize vapour flow
var Vap = V[1];
minimize VapourFlow: Vap;