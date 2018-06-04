#model of column with constant molar overflow and constant relative volatilites

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

var Vs >= 0; #vapor mass flow between stages
var y {S, C} >= 0, <= 1; #vapor concentrations between stages
subject to ClosureY {n in S}: 1 = sum{c in C} y[n,c];

var Ls >= 0; #liquid mass flow between stages
var x {S, C} >= 0, <= 1; #liquid concentrations between stages
subject to ClosureX {n in S}: 1 = sum{c in C} x[n,c]; 

var PC >= 0, <= Fall; #product mass flow after condenser
#same composition as y[NS,i], liquid
var PR >= 0, <= Fall; #product mass flow before reboiler
#same composition as x[1,i], liquid

param alpha {C}; #relative volatilities

subject to ComponentBalanceTopStage {c in C}:
Vs*y[NS-1,c] = PC*y[NS,c] + (Ls - Fall)*x[NS,c];
subject to ComponentBalance {n in Sm,c in C}:
Vs*y[n-1,c] + (Ls - sum{m in 2..n} F[m])*x[n+1,c] + F[n]*xF[c] = Vs*y[n,c] + (Ls - sum{m in 2..n-1} F[m])*x[n,c];
subject to ComponentBalanceBottomStage{c in C}:
Ls*x[2,c] = Vs*y[1,c] + PR*x[1,c];

subject to PositiveLiquid:
Ls - Fall >= 0;

subject to Thermodynamics {n in S, c in C}:
y[n,c]*(sum{d in C} alpha[d]*x[n,d]) =  alpha[c]*x[n,c];

#minimize vapour flow
var Vap = Vs;
minimize VapourFlow: Vap;





