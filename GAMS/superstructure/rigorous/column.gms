*model of one distillation column with n equilibrium stages

*choose solver
option NLP = KNITRO;

*choose mixture
$gdxin properties/BenzeneTolueneEthylbenzeneStyreneMethylstyrene.gdx

*choose subset of components;
set component /1, 2, 3/;
alias(component,c,d,e);

*choose column specific properties
set n set of stages /1*10/;
alias(n,m);
scalar NS number of stages;
NS = card(n);
scalar fn feed stage;
fn = 5;

*choose feed specific properties
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

Parameter Tlow lower boundary for temperature /100/;
Parameter Tup upper boundary for temperature /500/;
Parameter Qup upper boundary for heat duty /100000/;

Positive variable TF temperature of feed;
TF.lo = Tlow;
TF.up = Tup;
*starting point for T
TF.l = Tlow + (Tup-Tlow)/2;

Positive variable T(n) temperature on stage;
T.lo(n) = Tlow;
T.up(n) = Tup;
*starting point for T
T.l(n) = Tlow + (Tup-Tlow)/2;

Positive variable QR,QC 'heat duty at reboiler/condenser';
QR.up = Qup;
QC.up = Qup;
*starting point for QR/QC
QR.l = Qup/2;
QC.l = Qup/2;

Positive variable V(n) vapor mass flow;
*starting point for V
V.l(n) = Fall/2;

Positive variable y(n,c) vapor concentrations;
Equation ClosureY (n);
ClosureY(n) .. sum(c,y(n,c)) =e= 1;
*starting point for y
y.l(n,c) = xF(c);

Positive variable L(n) liquid mass flow;
*starting point for L
L.l(n) = Fall/2;

Positive variable x(n,c) liquid concentrations;
Equation ClosureX(n);
ClosureX(n) .. sum(c,x(n,c)) =e= 1;
*starting point for x
x.l(n,c) = xF(c);

Positive variable PR,PC 'product mass flow at reboiler/condenser';

Equation ComponentBalanceTopStage(c);
ComponentBalanceTopStage(c) ..
sum(n$(ord(n)=NS),V(n-1)*y(n-1,c)) =e=
sum(n$(ord(n)=NS),
L(n)*x(n,c) + PC*y(n,c));

Equation ComponentBalance(n,c);
ComponentBalance(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
V(n-1)*y(n-1,c) + L(n+1)*x(n+1,c) + F(n)*xF(c) =e=
V(n)*y(n,c) + L(n)*x(n,c);

Equation ComponentBalanceBottomStage(c);
ComponentBalanceBottomStage(c) ..
L('2')*x('2',c) =e=
V('1')*y('1',c) + PR*x('1',c);

Variable Href(c);
Equation DefineHref(c);
DefineHref(c) .. Href(c) =e= cH('1',c)*Tref + cH('2',c)*cH('3',c)/tanh(cH('3',c)/Tref) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/Tref);

Variable H(n,c);
Equation DefineH(n,c);
DefineH(n,c) .. H(n,c) =e= (cH('1',c)*T(n) + cH('2',c)*cH('3',c)/tanh(cH('3',c)/T(n)) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/T(n))-Href(c))/3.6e6;
Variable HF(c);
Equation DefineHF(c);
DefineHF(c) .. HF(c) =e= (cH('1',c)*TF + cH('2',c)*cH('3',c)/tanh(cH('3',c)/TF) + cH('4',c)*cH('5',c)*tanh(cH('5',c)/TF)-Href(c))/3.6e6;

Variable hL(n,c);
Equation DefinehL(n,c);
DefinehL(n,c) .. hL(n,c) =e= H(n,c) - (cHVL('1',c)*(1-T(n)/cHVL('2',c))**(cHVL('3',c) + T(n)/cHVL('2',c)*(cHVL('4',c) + T(n)/cHVL('2',c)* cHVL('5',c))))/3.6e6;
Variable hLF(c);
Equation DefinehLF(c);
DefinehLF(c) .. hLF(c) =e= HF(c) - (cHVL('1',c)*(1-TF/cHVL('2',c))**(cHVL('3',c) + TF/cHVL('2',c)*(cHVL('4',c) + TF/cHVL('2',c)* cHVL('5',c))))/3.6e6;

Equation EnergyBalanceTopStage;
EnergyBalanceTopStage ..
sum(c,sum(n$(ord(n)=NS),V(n-1)*y(n-1,c)*H(n-1,c))) =e=
sum(c,sum(n$(ord(n)=NS),L(n)*x(n,c)*hL(n,c) + PC*y(n,c)*hL(n,c)))+QC;

Equation EnergyBalance(n);
EnergyBalance(n)$(ord(n)>=2 and ord(n)<=NS-1) ..
sum(c,V(n-1)*y(n-1,c)*H(n-1,c) + L(n+1)*x(n+1,c)*hL(n+1,c) + F(n)*xF(c)*hLF(c)) =e=
sum(c,V(n)*y(n,c)*H(n,c) + L(n)*x(n,c)*hL(n,c));

Equation EnergyBalanceBottomStage;
EnergyBalanceBottomStage ..
sum(c,L('2')*x('2',c)*hL('2',c)) +QR =e=
sum(c,V('1')*y('1',c)*H('1',c) + PR*x('1',c)*hL('1',c));


Equation Thermodynamics(n,c);
Thermodynamics(n,c) .. P*y(n,c) =e=
x(n,c)*
exp(cPr('1',c) + cPr('2',c)/T(n) + cPr('3',c)*log(T(n)) + cPr('4',c)*T(n)**cPr('5',c))/cPr('6',c)*
exp(
sum(d,x(n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(n)))*(cGamma('2',d,c) + cGamma('3',d, c)/T(n)))/
sum(d,x(n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(n))))+
sum(d,
x(n,d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T(n)))*
(cGamma('2',c,d) + cGamma('3',c, d)/T(n) -
sum(e,x(n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(n)))*(cGamma('2',e,d) + cGamma('3',e,d)/T(n)))/
sum(e,x(n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(n)))))/
sum(e,x(n,e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T(n))))));

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


Equation LK;
LK ..  PC*sum(n$(ord(n)=NS),y(n,'1')) =g= 0.95*Fall*xF('1')
Equation HK;
HK .. PR*x('1','2') =g= 0.95*Fall*xF('2');

Variable Vap;
Equation DefineVap;
DefineVap .. Vap =e= V('1');


Model column / all /;

solve column using NLP minimizing Vap;




