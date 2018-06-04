$ontext
topology as variable
reformulation of discontinous constraints (max) with additional binary varibales
*too many discrete varibales for DEMO version, need to use platform
$offtext

$gdxin examples/Agrawal

Set component indices of components;
Set s size of power set of indices of free matrix entries;
*maximal value: 2^((n-2)*(n+1)/2)
$ load component,s

Scalar NC number of components;
NC = card(component);

alias(component,c,i,j,k,l,m,o);

Set U(i,j) indices of upper triangular matrix = possible mixtures;
U(i,j)$(ord(i) <= ord(j)) = yes;
Set hE(i,j,k,l) indices of horizontal edges = possible top product streams;
hE(i,j,k,l)$(U(i,j) and U(k,l) and ord(k) = ord(i) and (ord(l) >= ord(j)+1)) = yes;
Set dE(i,j,k,l) indices of diagonal edges = possible bottom product streams;
dE(i,j,k,l)$(U(i,j) and U(k,l) and (ord(k) = ord(i) + ord(l) - ord(j)) and (ord(l) >= ord(j)+1)) = yes;

Binary Variable x(i,j);
x.fx(i,j)$(ord(i) <= ord(j) and (ord(j) = 1 or ord(j) = NC)) = 1;
Binary Variable y(i,j,k,l);
Binary Variable z(i,j,k,l);

Equation LinkOuth(i,j);
LinkOuth(i,j)$(U(i,j) and ord(j) <= NC-1) ..
sum(l$hE(i,j,i,l),y(i,j,i,l)) =e= x(i,j);
*if entry of matrix exists if there is exactly one outgoing horizontal edge

Equation LinkOutd(i,j);
LinkOutd(i,j)$(U(i,j) and ord(j) <= NC-1) ..
sum((k,l)$dE(i,j,k,l),z(i,j,k,l)) =e= x(i,j);
* if entry of matrix exists if there is exactly one outgoing diagonal edge

Equation LinkInh(i,j);
LinkInh(i,j)$(U(i,j) and ord(j) >= 2)..
sum(l$hE(i,l,i,j),y(i,l,i,j)) =l= x(i,j);
*if there is an incoming horizontal edge, entry of matrix exists, number of incoming horizontal edges at most one

Equation LinkInd(i,j);
LinkInd(i,j)$(U(i,j) and ord(j) >= 2) ..
sum((k,l)$dE(k,l,i,j),z(k,l,i,j)) =l= x(i,j);
*if there is an incoming diagonal edge, entry of matrix exists, number of incoming diagonal edges at most one

Equation LinkIn(i,j);
LinkIn(i,j)$(U(i,j) and ord(j) >= 2) ..
sum(l$hE(i,l,i,j),y(i,l,i,j)) + sum((k,l)$dE(k,l,i,j),z(k,l,i,j)) =g= x(i,j);
*if entry of matrix exists, there is (at least) one incoming edge

Equation NondecreasingNumbers(i,j);
NondecreasingNumbers(i,j)$(U(i,j) and ord(j) <= NC-1) ..
sum(l$hE(i,j,i,l),y(i,j,i,l)*(NC+1-ord(l))) + sum((k,l)$dE(i,j,k,l),z(i,j,k,l)*(NC+1-ord(l))) =g= x(i,j)*(NC+1-ord(j));
*the numbers of components does not decrease from one entry to its sucessors

Parameter F feed flow;
Parameter xF(component) feed concentrations;
Parameter alpha(component) relative volatilites;

$load F, xF, alpha

Variable LK(i,j) light key komponent;
Equation DefineLK(i,j);
DefineLK(i,j)$(U(i,j) and ord(j) <= NC-1) ..
LK(i,j) =e= ord(i) + sum((k,l)$dE(i,j,k,l),z(i,j,k,l)*(ord(l)-ord(j))) - 1;
Variable HK(i,j);
Equation DefineHK(i,j) heavy key component;
DefineHK(i,j)$(U(i,j) and ord(j) <= NC-1) ..
HK(i,j) =e= ord(i) + NC - ord(j) - sum(l$hE(i,j,i,l),y(i,j,i,l)*(ord(l)-ord(j))) + 1;

Positive variable V(i,j) vapor flow in top section of column  ;
Equation ColumnSectionExists(i,j);
ColumnSectionExists(i,j)$(U(i,j) and ord(j)<= NC-1) ..
V(i,j) =l= 100*F*x(i,j);
*if column does not exist, vapor flow is zero

Positive variable D(i,j,k,l,c) componentwise top product stream;
Equation TopExists(i,j,k,l);
TopExists(i,j,k,l)$hE(i,j,k,l) ..
sum(c,D(i,j,k,l,c)) =l= F*y(i,j,k,l);
*if product stream does no exist, componentwise zero

Positive variable B(i,j,k,l,component) componentwise bottom product stream;
Equation BottomExists(i,j,k,l);
BottomExists(i,j,k,l)$dE(i,j,k,l) ..
sum(c,B(i,j,k,l,c)) =l= F*z(i,j,k,l);
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

Positive variable x1p(i,j,c), x1n(i,j,c);
Binary variable x1b(i,j,c);
*1 if c <= LK(i,j)-1, else 0
Equation x1e1(i,j,component), x1e2(i,j,c), x1e3(i,j,c);
x1e1(i,j,c) .. x1p(i,j,c) + x1n(i,j,c) =e= LK(i,j) - ord(c) + 1;
x1e2(i,j,c) .. x1p(i,j,c) =l= NC * x1b(i,j,c);
x1e3(i,j,c) .. x1n(i,j,c) =l= NC * (1-x1b(i,j,c));


Equation noLKinBottom(i,j,c);
noLKinBottom(i,j,c)$(U(i,j) and ord(j) <= NC-1) ..
x1b(i,j,c)*(sum((k,l)$dE(i,j,k,l),B(i,j,i,l,c))) =e= 0;
*all components lighter or equal light key do not go out in bottom product

Positive variable x2p(i,j,c), x2n(i,j,c);
Binary variable x2b(i,j,c);
Equation x2e1(i,j,c), x2e2(i,j,c), x2e3(i,j,c);
x2e1(i,j,c) .. x2p(i,j,c) + x2n(i,j,c) =e= ord(c) - HK(i,j);
x2e2(i,j,c) .. x2p(i,j,c) =l= NC * x2b(i,j,c);
x2e3(i,j,c) .. x2n(i,j,c) =l= NC * (1-x2b(i,j,c));

Equation noHKinTop(i,j,c);
noHKinTop(i,j,c)$(U(i,j) and ord(j) <= NC-1) ..
x2b(i,j,c)*(sum(l$hE(i,j,i,l),D(i,j,i,l,c))) =e= 0;
*all components heavier or equal heavy key do not go out in top product

Variable theta(i,j,c);
* Underwood roots in column
Equation ThetaLowerBound(i,j,c);
ThetaLowerBound(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) >= ord(i) and (ord(c) <= NC + ord(i) - ord(j) -1)) ..
theta(i,j,c) =g= alpha(c + 1);
Equation ThetaUpperBound (i,j,c);
ThetaUpperBound(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) >= ord(i) and (ord(c) <= NC + ord(i) - ord(j) -1)) ..
theta(i,j,c) =l= alpha(c);

Equation FirstUnderwoodFeed(c);
FirstUnderwoodFeed(c)$(ord(c) <= NC-1) ..
sum(m, alpha(m)*F*xF(m)/(alpha(m)-theta('1','1',c))) =e= 0;

Equation FirstUnderwood(i,j,c);
FirstUnderwood(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) >= ord(i) and (ord(c) <= NC + ord(i) - ord(j) -1)) ..
x(i,j)*(sum(m,alpha(m)*(sum(l$hE(i,l,i,j),D(i,l,i,j,m)) + sum((k,l)$dE(k,l,i,j),B(k,l,i,j,m)))/(alpha(m)-theta(i,j,c)))) =e= 0

Equation SecondUnderwood(i,j,c);
SecondUnderwood(i,j,c)$(U(i,j) and ord(j) <= NC-1 and ord(c) >= ord(i) and (ord(c) <= NC + ord(i) - ord(j) -1)) ..
x(i,j)*x1b(i,j,c)*x2b(i,j,c)*(V(i,j)-sum(m,alpha(m)*sum(l$hE(i,j,i,l),D(i,j,i,l,m))/(alpha(m)-theta(i,j,c)))) =g= 0;

Equation Equilibrium(i,j,k,l,m,c);
Equilibrium(i,j,k,l,m,c)$(U(i,j) and ord(j) >= 2 and ord(j) <= NC-1 and hE(i,l,i,j) and dE(k,m,i,j)) ..
x(i,j)*y(i,l,i,j)*z(k,m,i,j)*(D(i,l,i,j,c)*sum(o,B(k,m,i,j,o)*alpha(o)) - B(k,m,i,j,c)*alpha(c)*sum(o,D(i,l,i,j,o))) =e= 0;
*if product is produced by two columns, it must be a sidestream -> vapor(top product of lower column) and liquid(bottom product of upper column) in equilibrium

Equation IntegratedColumn(i,j,k,l,m);
IntegratedColumn(i,j,k,l,m)$(U(i,j) and ord(j) >= 2 and hE(i,l,i,j) and dE(k,m,i,j)) ..
x(i,j)*y(i,l,i,j)*z(k,m,i,j)*(V(i,l) - V(k,m)) =e= 0;
*if product is produced by two columns, they must actually be one column -> vapor flows in both columns are equal

Variable Reboiler;
Equation DefineReboiler;
DefineReboiler .. Reboiler =e=
sum((i,j)$(U(i,j) and ord(j) <= NC-1), V(i,j)) -
sum((i,j)$U(i,j), sum(l$hE(i,l,i,j), sum((k,m)$dE(k,m,i,j),x(i,j)*y(i,l,i,j)*z(k,m,i,j)*V(k,m))));
*if product is produced by two columns, only count the lowest reboiler duty

Model hybridmethod /all/;

Solve hybridmethod using MINLP minimizing Reboiler;

