$ontext
for fixed topology, given by x
calculate Underwood's method as NLP
$offtext

$gdxin examples/Agrawal

Set component indices of components;
$ load component

Scalar NC number of components;
NC = card(component);

alias(component,c,i,j,k,l,m,o);

Set U(i,j) indices of upper triangular matrix = possible mixtures;
U(i,j)$(ord(i) <= ord(j)) = yes;
Set hE(i,j,k,l) indices of horizontal edges = possible top product streams;
hE(i,j,k,l)$(U(i,j) and U(k,l) and ord(k) = ord(i) and (ord(l) >= ord(j)+1)) = yes;
Set dE(i,j,k,l) indices of diagonal edges = possible bottom product streams;
dE(i,j,k,l)$(U(i,j) and U(k,l) and (ord(k) = ord(i) + ord(l) - ord(j)) and (ord(l) >= ord(j)+1)) = yes;

*fix topology manually
Parameter x(i,j) / 1.1 1, 1.2 1, 1.3 1, 2.2 1, 2.3 1, 3.3 1/;
*calculate y and z from fixed x
Parameter y(i,j,k,l);
y(i,j,i,l)$(x(i,j) = 1 and hE(i,j,i,l) and x(i,l)=1 and ord(l)=smin(o$(hE(i,j,i,o) and x(i,o)=1),ord(o))) = 1;
Parameter z(i,j,k,l);
z(i,j,k,l)$(x(i,j) = 1 and dE(i,j,k,l) and x(k,l)=1 and ord(l)=smin((m,o)$(dE(i,j,m,o) and x(m,o)=1),ord(o))) = 1;

Parameter LK(i,j) light key component;
LK(i,j)$(U(i,j) and ord(j) <= NC-1) = ord(i) + sum((k,l)$dE(i,j,k,l),z(i,j,k,l)*(ord(l)-ord(j))) - 1;
Parameter HK(i,j) heavy key component;
HK(i,j)$(U(i,j) and ord(j) <= NC-1) = ord(i) + NC - ord(j) - sum(l$hE(i,j,i,l),y(i,j,i,l)*(ord(l)-ord(j))) + 1;

Parameter F feed flow;
Parameter xF(c) feed concentrations;
Parameter alpha(c) relative volatilites;

$load F, xF, alpha

Positive variable V(i,j) vapor flow in top section of column;
Equation ColumnSectionExists(i,j);
ColumnSectionExists(i,j)$(U(i,j) and (x(i,j)=0 or ord(j)= NC)) ..
V(i,j) =e= 0;
*if column does not exist, vapor flow is zero

Positive variable D(i,j,k,l,c) componentwise top product stream;
Equation TopExists(i,j,k,l);
TopExists(i,j,k,l)$(hE(i,j,k,l) and y(i,j,k,l)=0) ..
sum(c,D(i,j,k,l,c)) =e= 0;
*if product stream does no exist, componentwise zero

Positive variable B(i,j,k,l,c) componentwise bottom product stream;
Equation BottomExists(i,j,k,l);
BottomExists(i,j,k,l)$(dE(i,j,k,l) and z(i,j,k,l)=0) ..
sum(c,B(i,j,k,l,c)) =e= 0;
*if product stream does no exist, componentwise zero

Equation ComponentBalanceFeed(c);
ComponentBalanceFeed(c) ..
F*xF(c) =e=
sum(l$hE('1','1','1',l),D('1','1','1',l,c)) + sum(l$dE('1','1',l,l),B('1','1',l,l,c));

Equation ComponentBalance(i,j,c);
ComponentBalance(i,j,c)$(U(i,j) and ord(j) >= 2 and ord(j) <= NC-1) ..
sum(l$hE(i,l,i,j),D(i,l,i,j,c)) + sum((k,l)$dE(k,l,i,j),B(k,l,i,j,c)) =e=
sum(l$hE(i,j,i,l),D(i,j,i,l,c)) + sum((k,l)$dE(i,j,k,l),B(i,j,k,l,c));
*all ingoing component flow must go out

Equation noLKinBottom(i,j,c);
noLKinBottom(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) <= LK('1','1')) ..
sum((l,k)$dE(i,j,k,l),B(i,j,k,l,c)) =e= 0;
*all components lighter or equal light key do not go out in bottom product

Equation noHKinTop(i,j,c);
noHKinTop(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) >= HK(i,j)) ..
sum(l$hE(i,j,i,l),D(i,j,i,l,c)) =e= 0;
*all components heavier or equal heavy key do not go out in top product

Variable theta(i,j,c);
* Underwood roots in column
Equation ThetaLowerBound(i,j,c);
ThetaLowerBound(i,j,c)$(U(i,j) and x(i,j)=1 and ord(j) <= NC-1 and ord(c) >= LK(i,j) and ord(c) <= HK(i,j)-1) ..
theta(i,j,c) =g= alpha(c + 1);
Equation ThetaUpperBound (i,j,c);
ThetaUpperBound(i,j,c)$(U(i,j) and x(i,j)=1 and ord(j) <= NC-1 and ord(c) >= LK(i,j) and ord(c) <= HK(i,j)-1) ..
theta(i,j,c) =l= alpha(c);

Equation FirstUnderwoodFeed(c);
FirstUnderwoodFeed(c)$(ord(c) >= LK('1','1') and ord(c) <= HK('1','1')-1) ..
sum(m, alpha(m)*F*xF(m)/(alpha(m)-theta('1','1',c))) =e= 0;
Equation FirstUnderwood(i,j,c);

FirstUnderwood(i,j,c)$(U(i,j) and x(i,j)=1 and ord(j) <= NC-1 and ord(c) >= LK(i,j) and ord(c) <= HK(i,j)-1) ..
sum(m,alpha(m)*(sum(l$hE(i,l,i,j),D(i,l,i,j,m)) + sum((k,l)$dE(k,l,i,j),B(k,l,i,j,m)))/(alpha(m)-theta(i,j,c))) =e= 0;

Equation SecondUnderwood(i,j,c);
SecondUnderwood(i,j,c)$(U(i,j) and x(i,j)=1 and ord(j) <= NC-1 and ord(c) >= LK(i,j) and ord(c) <= HK(i,j)-1) ..
V(i,j)-sum(m,alpha(m)*sum(l$hE(i,j,i,l),D(i,j,i,l,m))/(alpha(m)-theta(i,j,c))) =g= 0;

Equation Equilibrium(i,j,k,l,m,c);
Equilibrium(i,j,k,l,m,c)$(U(i,j) and x(i,j)=1 and ord(j) >= 2 and ord(j) <= NC-1 and hE(i,l,i,j) and dE(k,m,i,j) and y(i,l,i,j)=1 and z(k,m,i,j)=1) ..
(D(i,l,i,j,c)*sum(o,B(k,m,i,j,o)*alpha(o)) - B(k,m,i,j,c)*alpha(c)*sum(o,D(i,l,i,j,o))) =e= 0;
*if product is produced by two columns, it must be a sidestream -> vapor(top product of lower column) and liquid(bottom product of upper column) in equilibrium

Equation IntegratedColumn(i,j,k,l,m);
IntegratedColumn(i,j,k,l,m)$(U(i,j) and x(i,j)=1 and ord(j) >= 2 and hE(i,l,i,j) and dE(k,m,i,j) and y(i,l,i,j) = 1 and z(k,m,i,j)=1) ..
V(i,l) - V(k,m) =e= 0;
*if product is produced by two columns, they must actually be one column -> vapor flows in both columns are equal

Variable Reboiler;
Equation DefineReboiler;
DefineReboiler .. Reboiler =e=
sum((i,j)$(U(i,j) and ord(j) <= NC-1), V(i,j)) -
sum((i,j)$U(i,j), sum(l$hE(i,l,i,j), sum((k,m)$dE(k,m,i,j),x(i,j)*y(i,l,i,j)*z(k,m,i,j)*V(k,m))));
*if product is produced by two columns, only count the vapour flow once

Model fixedTopology /all/;

Option NLP  = KNITRO;

solve fixedTopology using NLP minimizing Reboiler;


