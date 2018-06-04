*model of one distillation column with n equilibrium stages

*choose solver
option NLP = KNITRO;

*choose mixture
$gdxin properties/BenzeneTolueneEthylbenzeneStyreneMethylstyrene.gdx

*choose subset of components;
set component /1, 2, 3/;
alias(component,c,d);

*choose column specific properties
set n set of stages /1*15/;
alias(n,m);
scalar NS number of stages;
NS = card(n);
scalar fn feed stage;
fn = 8;

*choose feed specific properties
Parameter Fall feed flow /100/;
Parameter xF (c) feed concentrations /1 0.3, 2 0.4, 3 0.3/;

Parameter alpha (c) relative volatilities;
$load alpha
$gdxin

Positive variable F(n) feed flow on each stage;
F.fx(n)$(ord(n)<>fn) = 0;
F.fx(n)$(ord(n)=fn) = Fall;

Positive variable Vs vapor mass flow stripping;

Positive variable y(n,c) vapor concentrations;
Equation closureY (n);
closureY(n) .. sum(c,y(n,c)) =e= 1;
y.l(n,c) = xF(c);

Positive variable Ls liquid mass flow rectifying;

Positive variable x(n,c) liquid concentrations;
Equation closureX(n);
closureX(n) .. sum(c,x(n,c)) =e= 1;
x.l(n,c) = xF(c);

Positive variable PR,PC 'product mass flow at reboiler/condenser';

Equation ComponentBalanceTopStage(c);
ComponentBalanceTopStage(c) ..
sum(n$(ord(n)=NS),Vs*y(n-1,c)) =e=
sum(n$(ord(n)=NS),
(Ls-Fall)*x(n,c) + PC*y(n,c));

Equation ComponentBalance(n,c);
ComponentBalance(n,c)$(ord(n)>=2 and ord(n)<=NS-1) ..
Vs*y(n-1,c) + (Ls-sum(m$(ord(m)<=ord(n)),F(m)))*x(n+1,c) + F(n)*xF(c) =e=
Vs*y(n,c) + (Ls-sum(m$(ord(m)<=ord(n)-1),F(m)))*x(n,c);

Equation ComponentBalanceBottomStage(c);
ComponentBalanceBottomStage(c) ..
Ls*x('2',c) =e=
Vs*y('1',c) + PR*x('1',c);

Equation PositiveLiquid;
PositiveLiquid .. Ls - Fall =g= 0;

Equation Equilibrium(n,c);
Equilibrium(n,c) ..
y(n,c) * sum(d,alpha(d)*x(n,d)) =e= alpha(c)*x(n,c);

Equation LK;
LK ..  PC*sum(n$(ord(n)=NS),y(n,'1')) =g= 0.95*Fall*xF('1')
Equation HK;
HK .. PR*x('1','2') =g= 0.95*Fall*xF('2');

Variable Vap;
Equation DefineVap;
DefineVap .. Vap =e= Vs;

Model column / all /;

solve column using NLP minimizing Vap;





