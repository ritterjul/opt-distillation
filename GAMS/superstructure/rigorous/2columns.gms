*model for a superstructure of 2 distillation columns with n equilibrium stages

*choose solver
option NLP = KNITRO;
option MINLP = KNITRO;
option MPEC = KNITRO;

*choose mixture
$gdxin properties/BenzeneTolueneEthylbenzeneStyreneMethylstyrene.gdx

*choose subset of components;
set component /1, 2, 3/;
alias(component,c,d,e);

*choose column specific properties
set i set of columns /1*2/;
set n set of stages /1*5/;
alias(n,m);
scalar NS number of stages;
NS = card(n);
scalar fn feed stage;
fn = 3;
scalar ccn connecting vapour stage;
ccn = 4;
scalar crn connecting liquid stage;
crn = 2;

*feed specific properties
Parameter Fall feed flow /100/;
Parameter xF (c) feed concentrations /1 0.3, 2 0.4, 3 0.3/;

*mixture specific model constants
Parameter P pressure;
Parameter Tref reference temperature;
set index /1*6/;
Parameter cPr (index,c) model constants for vapor pressure;
set index2 /1*5/;
Parameter cH (index2, c) model constants for enthalpy of vapour;
Parameter cHVL (index2, c) model constants for heat of vaporisation;
set index3 /1*3/;
Parameter cGamma (index3, c, d) model constants for activity coefficient;
$load P, Tref, cPr, cH, cHVL, cGamma
$gdxin

Positive variable F(n) feed flow on each stage;
F.fx(n)$(ord(n)<>fn) = 0;
F.fx(n)$(ord(n)=fn) = Fall;

*column specific properties
Parameter Tlow lower boundary for temperature /275/;
Parameter Tup upper boundary for temperature /500/;
Parameter Qup upper boundary for heat duty /100000/;

Positive variable T(i,n) temperature;
T.lo(i,n) = Tlow;
T.up(i,n) = Tup;
*starting point for T
T.l(i,n) = Tlow + (Tup-Tlow)/2;
Positive variable TF temperature of feed;
*starting point for TF
TF.l = Tlow + (Tup-Tlow)/2;

Positive variable QR(i),QC(i) 'heat duty at reboiler/condenser';
QR.up(i) = Qup;
QC.up(i) = Qup;
*starting point for QR/QC
QR.l(i) = Qup/2;
QC.l(i) = Qup/2;

Positive variable V(i,n) vapor mass flow;
*starting point for V
V.l(i,n) = Fall;

Positive variable y (i,n,c) vapor concentrations;
Equation closureY (i,n);
closureY(i,n) .. sum(c,y(i,n,c)) =e= 1;
*starting point for y
y.l(i,n,c) = xF(c);

Positive variable L (i,n) liquid mass flow;
*starting point for L
L.l(i,n) = Fall;

Positive variable x(i,n,c) liquid concentrations;
Equation closureX (i,n);
closureX(i,n) .. sum(c,x(i,n,c)) =e= 1;
*starting point for x
x.l(i,n,c) = xF(c);

Positive variable PR1,PR2,PC1,PC2 'product mass flow at reboiler/condenser';
Positive variable PS(n) 'product mass flow at side draw';
Positive variable CR(n),CC(n) 'connecting mass flow at reboiler/condenser';

PS.fx(n)$(ord(n)<>fn) = 0;
CR.fx(n)$(ord(n)<>crn) = 0;
CC.fx(n)$(ord(n)<>ccn) = 0;

Equation ComponentBalanceTopStage1(c);
ComponentBalanceTopStage1(c) ..
sum(n$(ord(n)=NS),V('1',n-1)*y('1',n-1,c)) =e=
sum(n$(ord(n)=NS),
L('1',n)*x('1',n,c) + (PC1+sum(m,CC(m)))*y('1',n,c));

Equation ComponentBalance1(n,c);
ComponentBalance1(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
V('1',n-1)*y('1',n-1,c) + L('1',n+1)*x('1',n+1,c) + F(n)*xF(c) =e=
V('1',n)*y('1',n,c) + L('1',n)*x('1',n,c);

Equation ComponentBalanceBottomStage1(c);
ComponentBalanceBottomStage1(c) ..
L('1','2')*x('1','2',c) =e=
V('1','1')*y('1','1',c) + (PR1+sum(m,CR(m)))*x('1','1',c);

Equation ComponentBalanceTopStage2(c);
ComponentBalanceTopStage2(c) ..
sum(n$(ord(n)=NS),V('2',n-1)*y('2',n-1,c)) =e=
sum(n$(ord(n)=NS),
L('2',n)*x('2',n,c) + PC2*y('2',n,c));

Equation ComponentBalance2(n,c);
ComponentBalance2(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
V('2',n-1)*y('2',n-1,c) + L('2',n+1)*x('2',n+1,c) + CC(n)*sum(m$(ord(m)=NS),y('1',m,c)) + CR(n)*x('1','1',c) =e=
V('2',n)*y('2',n,c) + L('2',n)*x('2',n,c) + PS(n)*x('2',n,c);

Equation ComponentBalanceBottomStage2(c);
ComponentBalanceBottomStage2(c) ..
L('2','2')*x('2','2',c) =e=
V('2','1')*y('2','1',c) + PR2*x('2','1',c);

Variable Href(c);
Equation DefineHref(c);
DefineHref(c) .. Href(c) =e= cH('1',c)*Tref + cH('2',c)*cH('3',c)/tanh(cH('3',c)/Tref) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/Tref);

Variable H(i,n,c);
Equation DefineH(i,n,c);
DefineH(i,n,c) .. H(i,n,c) =e= (cH('1',c)*T(i,n) + cH('2',c)*cH('3',c)/tanh(cH('3',c)/T(i,n)) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/T(i,n))-Href(c))/3.6e6;
Variable HF(c);
Equation DefineHF(c);
DefineHF(c) .. HF(c) =e= (cH('1',c)*TF + cH('2',c)*cH('3',c)/tanh(cH('3',c)/TF) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/TF)-Href(c))/3.6e6;

Variable hL(i,n,c);
Equation DefinehL(i,n,c);
DefinehL(i,n,c) .. hL(i,n,c) =e= H(i,n,c) - (cHVL('1',c)*(1-T(i,n)/cHVL('2',c))**(cHVL('3',c) + T(i,n)/cHVL('2',c)*(cHVL('4',c) + T(i,n)/cHVL('2',c)* cHVL('5',c))))/3.6e6;
Variable hLF(c);
Equation DefinehLF(c);
DefinehLF(c) .. hLF(c) =e= HF(c) - (cHVL('1',c)*(1-TF/cHVL('2',c))**(cHVL('3',c) + TF/cHVL('2',c)*(cHVL('4',c) + TF/cHVL('2',c)* cHVL('5',c))))/3.6e6;

Equation EnergyBalanceTopStage1;
EnergyBalanceTopStage1 ..
sum(c,sum(n$(ord(n)=NS),V('1',n-1)*y('1',n-1,c)*H('1',n-1,c))) =e=
sum(c,sum(n$(ord(n)=NS),L('1',n)*x('1',n,c)*hL('1',n,c) + (PC1+sum(m,CC(m)))*y('1',n,c)*hL('1',n,c)))+QC('1');

Equation EnergyBalance1(n);
EnergyBalance1(n)$(ord(n)>=2 and ord(n)<=NS-1) ..
sum(c,V('1',n-1)*y('1',n-1,c)*H('1',n-1,c) + L('1',n+1)*x('1',n+1,c)*hL('1',n+1,c) + F(n)*xF(c)*hLF(c)) =e=
sum(c,V('1',n)*y('1',n,c)*H('1',n,c) + L('1',n)*x('1',n,c)*hL('1',n,c));

Equation EnergyBalanceBottomStage1;
EnergyBalanceBottomStage1 ..
sum(c,L('1','2')*x('1','2',c)*hL('1','2',c)) +QR('1') =e=
sum(c,V('1','1')*y('1','1',c)*H('1','1',c) + (PR1+sum(m,CR(m)))*x('1','1',c)*hL('1','1',c));

Equation EnergyBalanceTopStage2;
EnergyBalanceTopStage2 ..
sum(c,sum(n$(ord(n)=NS),V('2',n-1)*y('2',n-1,c)*H('2',n-1,c))) =e=
sum(c,sum(n$(ord(n)=NS),L('2',n)*x('2',n,c)*hL('2',n,c) + PC2*y('2',n,c)*hL('2',n,c)))+QC('2');

Equation EnergyBalance2(n);
EnergyBalance2(n)$(ord(n)>=2 and ord(n)<=NS-1) ..
sum(c,V('2',n-1)*y('2',n-1,c)*H('2',n-1,c) + L('2',n+1)*x('2',n+1,c)*hL('2',n+1,c) + CC(n)*sum(m$(ord(m)=NS),y('1',m,c)*hL('1',m,c)) + CR(n)*x('1','1',c)*hL('1','1',c)) =e=
sum(c,V('2',n)*y('2',n,c)*H('2',n,c) + (L('2',n)+PS(n))*x('2',n,c)*hL('2',n,c));

Equation EnergyBalanceBottomStage2;
EnergyBalanceBottomStage2 ..
sum(c,L('2','2')*x('2','2',c)*hL('2','2',c)) +QR('2') =e=
sum(c,V('2','1')*y('2','1',c)*H('2','1',c) + PR2*x('2','1',c)*hL('2','1',c));

Equation Thermodynamics(i,n,c);
Thermodynamics(i,n,c) .. P*y(i,n,c) =e=
x(i,n,c)*
exp(cPr('1',c) + cPr('2',c)/T(i,n) + cPr('3',c)*log(T(i,n)) + cPr('4',c)*T(i,n)**cPr('5',c))/cPr('6',c)*
exp(
sum(d,x(i,n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(i,n)))*(cGamma('2',d,c) + cGamma('3',d, c)/T(i,n)))/
sum(d,x(i,n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(i,n))))+
sum(d,
x(i,n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(i,n)))*
(cGamma('2',c,d) + cGamma('3',c, d)/T(i,n) -
sum(e,x(i,n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(i,n)))*(cGamma('2',e,d) + cGamma('3',e,d)/T(i,n)))/
sum(e,x(i,n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(i,n)))))/
sum(e,x(i,n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(i,n))))));
*x*PrS*gamma

Equation BoilingLiquidFeed;
BoilingLiquidFeed ..  P =e= sum(c,xF(c)*
exp(cPr('1',c) + cPr('2',c)/TF + cPr('3',c)*log(TF) + cPr('4',c)*TF**cPr('5',c))/cPr('6',c)*
exp(
sum(d,xF(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/TF))*(cGamma('2',d,c) + cGamma('3',d, c)/TF))/
sum(d,xF(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/TF)))+
sum(d,
xF(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/TF))*
(cGamma('2',c,d) + cGamma('3',c, d)/TF -
sum(e,xF(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/TF))*(cGamma('2',e,d) + cGamma('3',e,d)/TF))/
sum(e,xF(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/TF))))/
sum(e,xF(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/TF))))));

Variable Vap;
Equation DefineVap;
DefineVap .. Vap =e= V('1','1') + V('2','1');

Variable Sep;
Equation DefineSep;
DefineSep .. Sep =e=
PC1*sum(c,sum(n$(ord(n)=NS),y('1',n,c)**2))
+ PR1*sum(c,x('1','1',c)**2)
+ PC2*sum(c,sum(n$(ord(n)=NS),y('2',n,c)**2))
+ PR2*sum(c,x('2','1',c)**2)
+ sum(n,PS(n)*sum(c,x('2',n,c)**2));

*constrain overall separation
Equation RestrictSep;
RestrictSep .. Sep =g= 0.6*Fall;

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

Binary variables YPC1,YPR1,YPS;
*GDP
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






