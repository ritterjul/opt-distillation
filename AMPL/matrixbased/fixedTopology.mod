#Underwood model for fixed topology (given by matrix entries x)

#description of topology (fixed)
param NC integer; #size of matrix
set U := {i in 1..NC, j in 1..NC : i <= j}; #indices of upper triangle of matrix
param X {U} binary; #variables for entries in matrix

set E := {i in 1..NC-1, j in 1..NC-1, l in 1..NC: i <= j and l > j}; # indices of edges
param Y {(i,j,l) in E} = X[i,j]*X[i,l]*max(0,1-sum{m in j+1..l-1} X[i,m]); 
param Z {(i,j,l) in E} = X[i,j]*X[l-(j-i),l]*max(0,1-sum{m in j+1..l-1} X[m-(j-i),m]); 

param LK {(i,j) in U: j <> NC} = i + sum{k in 1..NC-j} Z[i,j,j+k] * k - 1;
param HK {(i,j) in U: j <> NC} = i + NC - j - sum{k in 1..NC-j} Y[i,j,j+k] * k + 1;

#description of feed and properties (fixed)
param Fall >= 0;
param xF {1..NC} >= 0, <= 1;
check: sum{m in 1..NC} xF[m] = 1;
param alpha {1..NC} >= 0;

var V {(i,j) in U: j <> NC} >= 0; # vapor flow in rectifying section
subject to ColumnSectionExists {(i,j) in U: j <> NC}: V[i,j] <= 100*Fall*X[i,j]; #if column does not exist, vapor flow is zero

var D {E, 1..NC} >= 0; #componentwise top product stream
subject to VaporFlowExists {(i,j,l) in E}: sum{m in 1..NC} D[i,j,l,m] <= Fall*Y[i,j,l]; #if stream does not exist, flow is zero
var B {E, 1..NC} >= 0; #componentwise bottom product stream 
subject to LiquidFlowExists {(i,j,l) in E}: sum{m in 1..NC} B[i,j,l,m] <= Fall*Z[i,j,l]; #if stream does not exist, flow is zero

#general component balances
subject to ComponentBalanceFeed {m in 1..NC}:
Fall*xF[m] = sum{l in 2..NC} D[1,1,l,m] + sum{l in 2..NC} B[1,1,l,m];
subject to ComponentBalance {j in 2..NC-1, i in 1..j, m in 1..NC}: 
sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m] = sum{l in j+1..NC} D[i,j,l,m] + sum{l in j+1..NC} B[i,j,l,m];
# all ingoing component flow must go out

#component balances given by separation
subject to noLKinBottom {j in 1..NC-1, i in 1..j, m in 1..LK[i,j]}: 
sum{l in j+1..NC} B[i,j,l,m] = 0;
#all components lighter or equal light key go out in top product
subject to noHKinTop {j in 1..NC-1, i in 1..j, m in HK[i,j]..NC}: 
sum{l in j+1..NC} D[i,j,l,m]= 0;
#all components heavier or equal heavy key do not go out in top product

#Underwood equations
var theta {i in 1..NC, j in 1..NC-1, r in 1..NC: i <= j and r >= i and r <= NC+i-j-1} >= alpha[r+1], <= alpha[r]; #Underwood roots in column

subject to FirstUnderwoodFeed {r in 1..NC-1}: 
sum{m in 1..NC} alpha[m]*Fall*xF[m]/(alpha[m]-theta[1,1,r]) = 0;
subject to FirstUnderwood {j in 2..NC-1, i in 1..j, r in i..NC+i-j-1: X[i,j]==1}: 
sum{m in 1..NC}(alpha[m]*(sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m]))/(alpha[m]-theta[i,j,r]) = 0; #all connecting streams liquid -> q=1

subject to SecondUnderwood {j in 1..NC-1, i in 1..j, r in LK[i,j]..HK[i,j]-1: X[i,j]==1}: 
V[i,j]-sum{m in 1..NC}(alpha[m]*(sum{l in j+1..NC} D[i,j,l,m]))/(alpha[m]-theta[i,j,r]) >= 0;

subject to IntegratedColumn {j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),l2,j] == 1}:
V[i,l1] = V[l2-(j-i),l2];  #if product is produced by two columns, they must actually be one column -> vapor flows in both columns are equal

subject to Equilibrium {j in 2..NC-1, i in 1..j, m in 1..NC, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),i,j] == 1}:  
D[i,l1,j,m]*(sum {o in 1..NC} B[l2-(j-i),l2,j,o]*alpha[o])- B[l2-(j-i),l2,j,m]*alpha[m]*(sum{o in 1..NC} D[i,l1,j,o]) = 0;
# if product is produced by two columns, it must be a sidestream -> vapor(top product of lower column) and liquid(bottom product of upper column) in equilibrium

minimize Vap: sum{j in 1..NC-1, i in 1..j} V[i,j] - sum{j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),l2,j] == 1} V[l2-(j-i),l2];
#if product is produced by two columns, only count the vapour in the lowest column section


/*
#in case of liquid and vapour connecting streams

redeclare subject to FirstUnderwood {j in 2..NC-1, i in 1..j, r in i..NC+i-j-1: X[i,j]=1}: 
sum{m in 1..NC}(alpha[m]*(sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m]))/(alpha[m]-theta[i,j,r]) = - sum{m in 1..NC}sum{l in i..j-1} D[i,l,j,m];

redeclare subject to IntegratedColumn {j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),l2,j] == 1}:
V[i,l1] = V[l2-(j-i),l2]-sum{m in 1..NC}sum{l in l2-(j-i)..l2-1} D[l2-(j-i),l,l2,m];

var Vs {i in 1..NC, j in 1..NC-1: i<=j} = V[i,j] - sum{m in 1..NC}sum{l in i..j-1} D[i,l,j,m];
#vapor flow in stripping section

redeclare minimize Vap: sum{j in 1..NC-1, i in 1..j} Vs[i,j] - sum{j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),l2,j] == 1} Vs[l2-(j-i),l2];
*/