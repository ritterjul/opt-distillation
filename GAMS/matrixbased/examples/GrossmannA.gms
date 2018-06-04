$ontext
example from Caballero, Grossmann: Design of distillation sequences:
from conventional to fullly thermally coupled distillation systems,
Computers & Chmeical Engineering
$offtext

Set component /1*5/;
Set s /1*512/;

Parameter F /100/;
Parameter xF(component) /1 0.04, 2 0.06, 3 0.5, 4 0.35, 5 0.05/;
Parameter alpha(component) /1 10.5, 2 4.04, 3 1.76, 4 1.31, 5 1.0/;

$gdxout examples/GrossmannA
$unload component, s, F, xF, alpha
$gdxout
