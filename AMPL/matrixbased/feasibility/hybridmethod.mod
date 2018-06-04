#description of the hybrid model with variables x,y,z and constraints

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

#redundant constraints
#subject to TotalCount: sum{j in 2..NC-1} sum{i in 1..j} X[i,j] >= NC-2; 
#configuration must involve at least NC-2 submixtures

#subject to PossiblePredecessor {j in 2..NC, i in 1..j}: sum{m in 1..j-i} X[i,j-m] + sum{m in 1..i-1} X[i-m,j-m] >= X[i,j]; 
#if entry of matrix exists, there must exist a possible predecessor