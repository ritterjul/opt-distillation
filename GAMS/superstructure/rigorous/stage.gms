*model of one equilibrium stage

*choose solver
option NLP = KNITRO;

*choose mixture
$gdxin properties/BenzeneTolueneEthylbenzeneStyreneMethylstyrene.gdx

*choose subset of components;
set component /1, 2, 3/;
alias(component,c,d,e);

*mixture specific model constants
Parameter P pressure;
set index /1*6/;
Parameter cPr (index,c) model constants for vapor pressure;
set index2 /1*5/;
Parameter cHVL (index, c) model constants for heat of vaporisation;
set index3 /1*3/;
Parameter cGamma (index3, c, d) model constants for activity coefficient;
$load P, cPr, cHVL, cGamma
$gdxin

*feed specific properties
Parameter F feed flow /10/;
Parameter xF (c) feed concentrations /1 0.3, 2 0.4, 3 0.3/;
* assume feed has same temperature as stage -> same enthalpy

*stage specific properties
Parameter Tlow lower boundary for temperature /100/;
Parameter Tup upper boundary for temperature /500/;
Parameter Qup upper boundary for heat duty /100/;

Positive variable T temperature;
T.lo = Tlow;
T.up = Tup;
*starting point for T
T.l = Tlow + (Tup-Tlow)/2;

Positive variable Q heat duty;
Q.up = Qup;
*starting point for Q
Q.l = Qup/2;

Positive variable V vapor mass flow;
*starting point for V
V.l = F/2;

Positive variable y(c) vapor concentrations;
Equation closureY;
closureY .. sum(c,y(c)) =e= 1;
*starting point for y
y.l(c) = xF(c);

Positive variable L liquid mass flow;
*starting point for L
L.l = F/2;

Positive variable x(c) liquid concentrations;
Equation closureX;
closureX .. sum(c,x(c)) =e= 1;
*starting point for x
x.l(c) = xF(c);

Equation ComponentBalance(c);
ComponentBalance(c) .. F*xF(c) =e=
V*y(c) + L*x(c);

Equation EnergyBalance;
EnergyBalance ..  Q =e=
sum(c,(F*xF(c)- L*x(c))*(cHVL('1',c)*(1-T/cHVL('2',c))**(cHVL('3',c) + T/cHVL('2',c)*(cHVL('4',c) + T/cHVL('2',c)* cHVL('5',c))))/3.6e6);
*(F*xF-L*x)*HVL

Equation Thermodynamics(c);
Thermodynamics(c) .. P*y(c) =e=
x(c)*
exp(cPr('1',c) + cPr('2',c)/T + cPr('3',c)*log(T) + cPr('4',c)*T**cPr('5',c))/cPr('6',c)*
exp(
sum(d,x(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T))*(cGamma('2',d,c) + cGamma('3',d, c)/T))/
sum(d,x(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T)))+
sum(d,
x(d)*exp(-cGamma('1',c,d)*(cGamma('2',c,d) + cGamma('3',c, d)/T))*
(cGamma('2',c,d) + cGamma('3',c, d)/T -
sum(e,x(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T))*(cGamma('2',e,d) + cGamma('3',e,d)/T))/
sum(e,x(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T))))/
sum(e,x(e)*exp(-cGamma('1',e,d)*(cGamma('2',e,d) + cGamma('3',e,d)/T)))));
*x*PrS*gamma

Equation EnoughVapor;
EnoughVapor .. V =g= 0.3*F;
Equation EnoughLiquid;
EnoughLiquid .. L =g= 0.3*F;

Equation ProductPurity;
ProductPurity ..  x('3') =g= 0.5;

Variable HeatDuty;
Equation DefineHeatDuty;
DefineHeatDuty .. HeatDuty =e= Q;

Model stage /all/;

solve stage using NLP minimizing HeatDuty;



