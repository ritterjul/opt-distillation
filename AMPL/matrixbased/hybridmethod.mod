#Underwood model for variable topology(given by  variables x,y,z) based on hybrid model

#description of topology
param NC integer; #size of matrix
set U := {i in 1..NC, j in 1..NC : i <= j}; #indices of upper triangle of matrix
var X {U} binary; #variables for entries in matrix
set E := {i in 1..NC-1, j in 1..NC-1, l in 1..NC: i <= j and l > j}; # indices of edges
var Y {E} binary; #variables for horizontal edges
var Z {E} binary; #variables for diagonal edges

subject to LinkInh {j in 2..NC, i in 1..j}: sum{l in i..j-1} Y[i,l,j] <= X[i,j]; 
#if there is an incoming horizontal edge, entry of matrix exists
subject to LinkInd {j in 2..NC, i in 1..j}: sum {l in j-i+1..j-1} Z[l-(j-i),l,j] <= X[i,j]; 
#if there is an incoming diagonal edge, entry of matrix exists
subject to LinkOuth {j in 1..NC-1, i in 1..j}: sum{l in j+1..NC} Y[i,j,l] = X[i,j]; 
#if entry of matrix exists if there is exactly one outgoing horizontal edge
subject to LinkOutd {j in 1..NC-1, i in 1..j}: sum{l in j+1..NC} Z[i,j,l] = X[i,j]; 
# if entry of matrix exists if there is exactly one outgoing diagonal edge
subject to LinkingIn {j in 2..NC, i in 1..j}: sum{l in i..j-1} Y[i,l,j] + sum {l in j-i+1..j-1} Z[l-(j-i),l,j] >= X[i,j]; 
#if entry of matrix exists, there is (at least) one incoming edge
subject to NondecreasingNumbers {j in 1..NC-2, i in 1..j}: sum{l in j+1..NC} (Y[i,j,l]+Z[i,j,l])*(NC+1-l) >= X[i,j]*(NC+1-j); 
#the NCumbers of components does NCot decrease from one entry to its sucessors

var LK {(i,j) in U: j <> NC} = i + sum{k in 1..NC-j} Z[i,j,j+k] * k - 1; #light key component
var HK {(i,j) in U: j <> NC} = i + NC - j - sum{k in 1..NC-j} Y[i,j,j+k] * k + 1; #heavy key component

#description of feed and properties (fixed)
param Fall >= 0;
param xF {1..NC} >= 0, <= 1;
check: sum{m in 1..NC} xF[m] = 1;
param alpha {1..NC} >= 0;

var V {(i,j) in U: j <> NC} >= 0; # vapor flow in rectifying section
subject to ColumnSectionExists {(i,j) in U: j <> NC}: V[i,j] <= 100*Fall*X[i,j]; #if column does NCot exist, vapor flow is zero

var D {E, 1..NC} >= 0; #componentwise top product stream
subject to VaporFlowExists {(i,j,l) in E}: sum{m in 1..NC} D[i,j,l,m] <= Fall*Y[i,j,l]; #if stream does NCot exist, flow is zero
var B {E, 1..NC} >= 0; #componentwise bottom product stream
subject to LiquidFlowExists {(i,j,l) in E}: sum{m in 1..NC} B[i,j,l,m] <= Fall*Z[i,j,l]; #if stream does NCot exist, flow is zero

#general component balances
subject to ComponentBalanceFeed {m in 1..NC}:
Fall*xF[m] = sum{l in 2..NC} D[1,1,l,m] + sum{l in 2..NC} B[1,1,l,m];
subject to ComponentBalance {j in 2..NC-1, i in 1..j, m in 1..NC}: 
sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m] = sum{l in j+1..NC} D[i,j,l,m] + sum{l in j+1..NC} B[i,j,l,m];
# all ingoing component flow must go out

#component balances given by separation
subject to noLKinBotom {j in 1..NC-1, i in 1..j, m in 1..NC}: 
max(LK[i,j]-m+1,0)*(sum{l in j+1..NC} B[i,j,l,m]) = 0;
#all components lighter or equal light key do not go out in bottom product
subject to noHKinTop {j in 1..NC-1, i in 1..j, m in 1..NC}: 
max(m-HK[i,j]+1,0)*(sum{l in j+1..NC} D[i,j,l,m]) = 0;
#all components heavier or equal heavy key do not go out in top product

var theta {i in 1..NC, j in 1..NC-1, r in 1..NC: i <= j and r >= i and r <= NC+i-j-1} >= alpha[r+1], <= alpha[r]; #Underwood roots in column

subject to FirstUnderwoodFeed {r in 1..NC-1}: 
sum{m in 1..NC} alpha[m]*Fall*xF[m]/(alpha[m]-theta[1,1,r]) = 0;
subject to FirstUnderwood {j in 2..NC-1, i in 1..j, r in i..NC+i-j-1}: 
X[i,j]*(sum{m in 1..NC}(alpha[m]*(sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m]))/(alpha[m]-theta[i,j,r])) = 0; #all connecting streams liquid -> q=1

subject to SecondUnderwood {j in 1..NC-1, i in 1..j, r in 1..NC: r >= i and r <= NC+i-j-1}: 
X[i,j]*max(r-LK[i,j]+1,HK[i,j]-r,0)*(V[i,j]-sum{m in 1..NC}(alpha[m]*sum{l in j+1..NC} D[i,j,l,m])/(alpha[m]-theta[i,j,r])) >= 0;

subject to IntegratedColumn {j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1}:
X[i,j]*Y[i,l1,j]*Z[l2-(j-i),l2,j]*(V[i,l1] - V[i+l2-j,l2]) = 0;  #if product is produced by two columns, they must actually be one column -> vapor flows in both columns are equal

subject to Equilibrium {j in 2..NC-1, i in 1..j, m in 1..NC, l1 in i..j-1, l2 in j-i+1..j-1}:  
X[i,j]*Y[i,l1,j]*Z[i+l2-j,l2,j]*(D[i,l1,j,m]*(sum {o in 1..NC} B[l2-(j-i),l2,j,o]*alpha[o]) - B[l2-(j-i),l2,j,m]*alpha[m]*(sum{o in 1..NC} D[i,l1,j,o])) = 0; 
#if product is produced by two columns, it must be a sidestream -> vapor(top product of lower column) and liquid(bottom product of upper column) in equilibrium

minimize Vap: sum{j in 1..NC-1, i in 1..j} V[i,j] - sum{j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1} X[i,j]*Y[i,l1,j]*Z[l2-(j-i),l2,j]*V[i+l2-j,l2];
#if product is produced by two columns, only count the lowest reboiler duty


/*
#in case of liquid and vapour connecting streams

redeclare subject to FirstUnderwood {j in 2..NC-1, i in 1..j, r in i..NC+i-j-1}: 
X[i,j]*(sum{m in 1..NC}(alpha[m]*(sum{l in i..j-1} D[i,l,j,m] + sum {l in j-i+1..j-1} B[l-(j-i),l,j,m]))/(alpha[m]-theta[i,j,r])-(sum{m in 1..NC}sum{l in i..j-1} D[i,l,j,m])) = 0;

redeclare subject to IntegratedColumn {j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1}:
X[i,j]*Y[i,l1,j]*Z[l2-(j-i),l2,j]*(V[i,l1]- V[l2-(j-i),l2]+sum{m in 1..NC}sum{l in i+l2-j..l2-1} D[i+l2-j,l,l2,m]) = 0;  

var Vs {i in 1..NC, j in 1..NC-1: i<=j} = V[i,j] - sum{m in 1..NC}sum{l in i..j-1} D[i,l,j,m];
#vapor flow in stripping section

redeclare minimize Vap: sum{j in 1..NC-1, i in 1..j} Vs[i,j] - sum{j in 2..NC, i in 1..j, l1 in i..j-1, l2 in j-i+1..j-1: X[i,j]== 1 and Y[i,l1,j] == 1 and  Z[l2-(j-i),l2,j] == 1} Vs[l2-(j-i),l2];
*/
