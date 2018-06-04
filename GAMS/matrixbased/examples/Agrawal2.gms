$ontext
example from Agrawal: Synthesis of Energy Efficient Distillation 
Configurations and New Dividing Wall Column, WCCE
$offtext

Set component /1*4/;
Set s /1*32/;

Parameter F /0.5/;
*NOT MENTIONED ON SLIDES
Parameter xF(component) /1 0.25, 2 0.25, 3 0.25, 4 0.25/;
Parameter alpha(component) /1 8.11, 2 3.53, 3 1.5, 4 1/;

$gdxout examples/Agrawal2
$unload component, s, F, xF, alpha
$gdxout
