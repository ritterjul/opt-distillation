$ontext
example from Caballero, Grossmann: Design of distillation sequences:
from conventional to fullly thermally coupled distillation systems,
Computers & Chmeical Engineering
$offtext

Set component /1*5/;
Set s /1*512/;

Parameter F /200/;
Parameter xF(component) /1 0.2, 2 0.3, 3 0.2, 4 0.2, 5 0.1/;
Parameter alpha(component) /1 4.1, 2 3.6, 3 2.1, 4 1.42, 5 1.0/;

$gdxout examples/GrossmannB
$unload component, s, F, xF, alpha
$gdxout
