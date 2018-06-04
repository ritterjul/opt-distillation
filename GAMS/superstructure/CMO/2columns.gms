*model for a superstructure of 2 distillation columns with n equilibrium stages

*choose solver
option NLP = KNITRO;

*choose mixture
$gdxin properties/BenzeneTolueneEthylbenzeneStyreneMethylstyrene.gdx

*choose subset of components;
set component /1, 2, 3/;
alias(component,c,d)

*choose column specific properties
set i set of columns /1*2/;
set n maximal set of stages /1*7/;
scalar NS;
NS = card(n);
alias(n,m);
scalar fn feed stage;
fn = 4;
scalar ccn connecting vapour stage;
ccn = 6;
scalar crn connecting liquid stage;
crn = 2;

*choose feed specific properties
Parameter Fall feed flow /100/;
Parameter xF (c) feed concentrations /1 0.3, 2 0.4, 3 0.3/;
Parameter alpha (c) relative volatilities;
$load alpha
$gdxin

Positive variable F(n) feed flow on each stage;
F.fx(n)$(ord(n)<>fn) = 0;
F.fx(n)$(ord(n)=fn) = Fall;

Positive variable Vs(i) vapor mass flow stripping;
Vs.l(i) = Fall;

Positive variable y (i,n,c) vapor concentrations;
Equation closureY (i,n);
closureY(i,n) .. sum(c,y(i,n,c)) =e= 1;
y.l(i,n,c) = xF(c);

Positive variable Ls(i) liquid mass flow stripping;
Ls.l(i) = Fall;

Positive variable x(i,n,c) liquid concentrations;
Equation closureX (i,n);
closureX(i,n) .. sum(c,x(i,n,c)) =e= 1;
x.l(i,n,c) = xF(c);

Positive variable PR1, PR2, PC1, PC2 'product mass flow at reboiler/condenser';
Positive variable PS(n) 'product mass flow at side draw'
Positive variable CR(n),CC(n) 'connecting mass flow at reboiler/condenser';

PS.fx(n)$(ord(n)<>fn) = 0;
CR.fx(n)$(ord(n)<>crn) = 0;
CC.fx(n)$(ord(n)<>ccn) = 0;

Equation ComponentBalanceTopStage1(c);
ComponentBalanceTopStage1(c) ..
sum(n$(ord(n)=NS),Vs('1')*y('1',n-1,c)) =e=
sum(n$(ord(n)=NS),(Ls('1')-Fall)*x('1',n,c) + (PC1 + sum(m,CC(m)))*y('1',n,c));

Equation ComponentBalance1(n,c);
ComponentBalance1(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
Vs('1')*y('1',n-1,c) + (Ls('1')-sum(m$(ord(m)<=ord(n)),F(m)))*x('1',n+1,c) =e=
Vs('1')*y('1',n,c) + (Ls('1')-sum(m$(ord(m)<=ord(n)-1),F(m)))*x('1',n,c);

Equation ComponentBalanceBottomStage1(c);
ComponentBalanceBottomStage1(c) ..
Ls('1')*x('1','2',c) =e=
Vs('1')*y('1','1',c) + (PR1 + sum(m,CR(m)))*x('1','1',c);

Equation ComponentBalanceTopStage2(c);
ComponentBalanceTopStage2(c) ..
sum(n$(ord(n)=NS),Vs('2')*y('2',n-1,c)) =e=
sum(n$(ord(n)=NS),PC2*y('2',n,c) + (Ls('2')-sum(m,CC(m)+CR(m)-PS(m)))*x('2',n,c));

Equation ComponentBalance2(n,c);
ComponentBalance2(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
Vs('2')*y('2',n-1,c) + (Ls('2') - sum(m$(ord(m)<=ord(n)),CC(m)+CR(m)-PS(m)))*x('2',n+1,c) + CC(n)*sum(m$(ord(m)=NS),y('1',m,c)) + CR(n)*x('1','1',c) =e=
Vs('2')*y('2',n,c) + (Ls('2') - sum(m$(ord(m)<=ord(n)),CC(m)+CR(m)-PS(m)))*x('2',n,c) + PS(n)*x('2',n,c);

Equation ComponentBalanceBottomStage2(c);
ComponentBalanceBottomStage2(c) ..
Ls('2')*x('2','2',c) =e=
Vs('2')*y('2','1',c) + PR2*x('2','1',c);

Equation PositiveLiquid1;
PositiveLiquid1 .. Ls('1') =g= Fall;

Equation PositiveLiquid2Strip;
PositiveLiquid2Strip .. Ls('2') =g= sum(m,CR(m));

Equation PositiveLiquid2Rec;
PositiveLiquid2Rec .. Ls('2') =g= sum(m,CC(m)+CR(m)-PS(m));

Equation Equilibrium(i,n,c);
Equilibrium(i,n,c) ..
y(i,n,c) * sum(d,alpha(d)*x(i,n,d)) =e= alpha(c)*x(i,n,c);

Variable Vap;
Equation DefineVap;
DefineVap .. Vap =e= Vs('1')+ Vs('2');

*get pure products
Variable Sep;
Equation DefineSep;
DefineSep .. Sep  =e=
PC1*sum(n$(ord(n)=NS),sum(c,y('1',n,c)**2))
+ PR1*sum(c,x('1','1',c)**2)
+ PC2*sum(n$(ord(n)=NS),sum(c,y('2',n,c)**2))
+ PR2*sum(c,x('2','1',c)**2)
+ sum(n,PS(n)*sum(c,x('2',n,c)**2));

*constrain overall separation
Equation RestrictSep;
RestrictSep .. Sep =g= 0.8*Fall;

Model columns / all /;
*solve columns using NLP minimizing Vap;

*restrict to 3 topologies
*MINLP
Binary variable yPC1,yPR1,yPS;
Equation Top, TopPC1, TopPR1,TopPS;
Top ..  1 =e= yPC1 + yPR1 + yPS;
TopPC1 .. 0 =g= (1-yPC1)*PC1 + yPC1*sum(n,CC(n));
TopPR1 .. 0 =g= (1-yPR1)*PR1 + yPR1*sum(n,CR(n));
TopPS .. 0 =g= (1-yPS)*sum(n,PS(n));

Model columnsMINLP /columns, Top, TopPC1,TopPR1,TopPS/;
*solve columnsMINLP using MINLP minimizing Vap;

*MPCC
Equation CCexists;
CCexists .. sum(n,CC(n)) =g= 0;
Equation CRexists;
CRexists .. sum(n,CR(n)) =g= 0;
Equation PSexists;
PSexists .. sum(n,PS(n)) =g= 0;
Positive Variable PCPR;
Equation DefinePCPR;
DefinePCPR .. PCPR =e= PC1 + PR1;
Equation PCPRnot;
PCPRnot .. PC1*PR1 =l= 0;

Model columnsMPCC /      columns,DefinePCPR, CCexists.PC1,CRexists.PR1,PSexists.PCPR,PCPRnot /;
*solve columnsMPCC using MPEC minimizing Vap;
*infeasible solution!

$ontext
*fix topology
PC1.fx = 0;
CR.fx(n) = 0;
PS.fx(n) = 0;

*solve columnsMPCC using MPEC minimizing Vap;
*works!
$offtext

*GDP 
Binary variables YPC1,YPR1,YPS;
Equation YPC1exists;
YPC1exists .. sum(n, CC(n)+PS(n))+PR1 =l= 0;
Equation YPR1exists;
YPR1exists .. sum(n, CR(n)+PS(n))+PC1 =l= 0;
Equation YPSexists;
YPSexists .. PR1+PC1 =l= 0;
Equation Dummy;
Dummy .. Sep =g= 0;

Logic equation LEq;
LEq .. (YPC1 and not YPR1 and not YPS) or  (not YPC1 and YPR1 and not YPS) or  (not YPC1 and not YPR1 and YPS);

Model columnsGDP /       columns, Dummy,
                         YPC1exists, YPR1exists, YPSexists,
                         LEq /;
File emp / '%emp.info%' /;
putclose emp
'disjunction YPC1 YPC1exists else Dummy' /
'disjunction YPR1 YPR1exists else Dummy' /
'disjunction YPS YPSexists else Dummy' /
;

solve columnsGDP using EMP minimizing Vap;
