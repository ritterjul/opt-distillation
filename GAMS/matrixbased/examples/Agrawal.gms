$ontext
example from Giridhar, Agrawal: Synthesis of distillation configurations II:
A search formulation for basic configurations, Computers and Chemical Engineering
$offtext

Set component /1*4/;
Set s /1*32/;

Parameter F /0.4/;
Parameter xF(component) /1 0.17757, 2 0.40458, 3 0.27708, 4 0.14077/;
Parameter alpha(component) /1 3.413, 2 2.146, 3 1.398, 4 1.0/;

$gdxout examples/Agrawal
$unload component, s, F, xF, alpha
$gdxout
