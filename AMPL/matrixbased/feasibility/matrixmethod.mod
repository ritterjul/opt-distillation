#description of the matrix model with variables x and constraints

param NC integer; #size of matrix

set U := {i in 1..NC, j in 1..NC : i <= j}; #indices of upper triangle of matrix

var X {U} binary; #variables for entries in matrix

subject to Predecessor {(i,j) in U: j in 2..NC}:
sum{m in 1..j-i} X[i,j-m] + sum{m in 1..i-1} X[i-m,j-m] >= X[i,j];

subject to NonDecreasingNumbers {(i,j) in U, k in 1..NC-j, l in 1..NC-j: j in 1..NC-2}:
X[i,j]*((k+1)*max(0,1-sum{m in 1..k} X[i,j+m]) + (l+1)*max(0,1-sum{m in 1..j} X[i+m,j+m]) + j) <= NC+1;