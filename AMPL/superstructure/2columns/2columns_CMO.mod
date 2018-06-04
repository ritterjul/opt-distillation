#model of superstructure for two columns with constant molar overflow and constant relative volatilities

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

var Vs {i in I} >= 0; #vapor mass flow rectifying section
var y {i in I, S[i], C} >= 0, <= 1; #vapor concentrations between stages
subject to ClosureY {i in I, n in S[i]}: 1 = sum{c in C} y[i,n,c];

var Ls {i in I} >= 0; #liquid mass flow stripping section
var x {i in I, S[i], C} >= 0, <= 1; #liquid concentrations between stages
subject to ClosureX {i in I, n in S[i]}: 1 = sum{c in C} x[i,n,c]; 

var CC {Sm[2]} >= 0, <= Fall; #connecting mass flow from before condenser of first column
#same composition as y[1,NS[1],i], vapour
var CR {Sm[2]} >= 0, <= Fall; #connecting mass flow from before reboiler of first column
#same composition as x[1,1,i], liquid
var PC {I} >= 0, <= Fall; #product mass flow after condenser
#same composition as y[i,NS[i],i], liquid
var PR {I} >= 0, <= Fall; #product mass flow before reboiler
#same composition as x[i,1,i], liquid
var PS {Sm[2]} >= 0, <= Fall; #side product mass flow at side
#same concentration as x[2,psn,i], liquid

param alpha {C}; #relative volatilities

subject to ComponentBalanceTopStage1 {c in C}:
Vs[1]*y[1,NS[1]-1,c] = (PC[1] + sum{m in Sm[2]} CC[m])*y[1,NS[1],c] + (Ls[1]-Fall)*x[1,NS[1],c];

subject to ComponentBalance1 {n in Sm[1],c in C}:
Vs[1]*y[1,n-1,c] + (Ls[1] - sum{m in 2..n} F[m])*x[1,n+1,c] + F[n]*xF[c] = Vs[1]*y[1,n,c] + (Ls[1] - sum{m in 2..n-1} F[m])*x[1,n,c];

subject to ComponentBalanceBottomStage1{c in C}:
Ls[1]*x[1,2,c] = Vs[1]*y[1,1,c] + (PR[1] + sum{m in Sm[2]} CR[m])*x[1,1,c];

subject to ComponentBalanceTopStage2{c in C}:
Vs[2]*y[2,NS[2]-1,c] = PC[2]*y[2,NS[2],c] + (Ls[2] - sum{m in Sm[2]} (CC[m]+CR[m]-PS[m]))*x[2,NS[2],c];

subject to ComponentBalance2 {n in Sm[2],c in C}:
Vs[2]*y[2,n-1,c] + (Ls[2] - sum{m in 2..n} (CC[m]+CR[m]-PS[m]))*x[2,n+1,c] + CC[n]*y[1,NS[1],c] + CR[n]*x[1,1,c] = Vs[2]*y[2,n,c] + (Ls[2] - sum{m in 2..n-1} (CC[m]+CR[m]-PS[m]))*x[2,n,c] + PS[n]*x[2,n,c];

subject to ComponentBalanceBottomStage2{c in C}:
Ls[2]*x[2,2,c] = Vs[2]*y[2,1,c] + PR[2]*x[2,1,c];

subject to IositiveLiquid1:
Ls[1] >= Fall;

subject to IositiveLiquid2Strip: 
Ls[2] >= sum{m in Sm[2]} CR[m];

subject to IositiveLiquid2Rec: 
Ls[2] >= sum{m in Sm[2]} (CC[m]+CR[m]-PS[m]);

subject to Thermodynamics {i in I, n in S[i], c in C}:
y[i,n,c]*(sum{d in C} alpha[d]*x[i,n,d]) =  alpha[c]*x[i,n,c];

var Vap = Vs[1] + Vs[2];
minimize VapourFlow: Vap;





