#description of the graph model with variables Y,Z and constraints

param NC integer; # siZe of matrix

set E := {i in 1..NC, j in 1..NC-1, l in 1..NC: i <= j and l > j}; # indices of edges

var Y {E} binary;
var Z {E} binary;

subject to hStart: sum{l in 2..NC} Y[1,1,l] = 1; #outgoing horizontal edge from start
subject to dStart: sum{l in 2..NC} Z[1,1,l] = 1; #outgoing diagonal edge from start

subject to hOuthIn {j in 2..NC-1, i in 1..j}: sum{l in i..j-1} Y[i,l,j] <= sum{l in j+1..NC} Y[i,j,l] ; 
#if there is an ingoing horiZontal edge, there is an outgoing horizontal edge
subject to dOuthIn {j in 2..NC-1, i in 1..j}: sum{l in i..j-1} Y[i,l,j] <= sum{l in j+1..NC} Z[i,j,l]; 
#if there is an ingoing horiZontal edge, there is an outgoing diagonal edge
subject to hOutdIn {j in 2..NC-1, i in 1..j}: sum{l in j-i+1..j-1} Z[l-(j-i),l,j] <= sum{l in j+1..NC} Y[i,j,l] ; 
#if there is an ingoing diagonal edge, there is an outgoing horizontal edge
subject to dOutdIn {j in 2..NC-1, i in 1..j}: sum{l in j-i+1..j-1} Z[l-(j-i),l,j] <= sum{l in j+1..NC} Z[i,j,l]; 
#if there is an ingoing diagonal edge, there is an outgoing diagonal edge

subject to InhOut {j in 2..NC-1, i in 1..j}: sum{l in i..j-1} Y[i,l,j] + sum{l in j-i+1..j-1} Z[l-(j-i),l,j] >= sum{l in j+1..NC} Y[i,j,l]; 
#if there is an outgoing horizontal edge, there is an ingoing edge
subject to IndOut {j in 2..NC-1, i in 1..j}: sum{l in i..j-1} Y[i,l,j] + sum{l in j-i+1..j-1} Z[l-(j-i),l,j] >= sum{l in j+1..NC} Z[i,j,l]; 
#if there is an outgoing diagonal edge, there is an ingoing edge

subject to End {i in 1..NC}: sum{l in i..NC-1} Y[i,l,NC] + sum{l in NC-i+1..NC-1} Z[l-(NC-i),l,NC] >= 1; 
#ingoing edge at end

subject to redunIn: sum{l in 2..NC} (Y[1,1,l]+Z[1,1,l])*(NC+1-l) >= NC; #all components go out at start

subject to redunh {j in 2..NC-1, i in 1..j}: sum{l in j+1..NC} (Y[i,j,l]+Z[i,j,l])*(NC+1-l) >= sum{l in i..j-1} Y[i,l,j]*(NC+1-j); 
#all components go out, if NCode is active
subject to redund {j in 2..NC-1, i in 1..j}: sum{l in j+1..NC} (Y[i,j,l]+Z[i,j,l])*(NC+1-l) >= sum{l in j-i+1..j-1} Z[l-(j-i),l,j]*(NC+1-j); 
#all components go out, if NCode is active

subject to basich {i in 1..NC}: sum{l in i..NC-1} Y[i,l,NC] <= 1; 
#there is at most one ingoing hoirzontal edge
subject to basicd {i in 1..NC}: sum{l in NC-i+1..NC-1} Z[l-(NC-i),l,NC] <= 1; 
#there is at most one ingoing diagonal edge