$ontext
example from Caballero, Grossmann: Design of distillation sequences:
from conventional to fullly thermally coupled distillation systems,
Computers & Chmeical Engineering
$offtext

Set component /1*5/;
Set s /1*512/;

Parameter F /200/;
Parameter xF(component) /1 0.1, 2 0.2, 3 0.3, 4 0.2, 5 0.2/;
Parameter alpha(component) /1 8.9, 2 5.7, 3 3.2, 4 1.55, 5 1.0/;

$gdxout examples/GrossmannC
$unload component, s, F, xF, alpha
$gdxout
