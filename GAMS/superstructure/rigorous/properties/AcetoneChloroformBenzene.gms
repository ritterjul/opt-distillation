*0- confidential
set component /1*3/;

Parameter P /0/;

set index /1*6/;

Table cPr(index,component)
                                 '2'     '1'     '3'
                        '1'      0       0       0
                        '2'      0       0       0
                        '3'      0       0       0
                        '4'      0       0       0
                        '5'      0       0       0
                        '6'      0       0       0
;

Parameter Tref /298/;

set index2 /1*5/;
Table cHV(index2,component)
                                 '2'     '1'     '3'
                        '1'      0       0       0
                        '2'      0       0       0
                        '3'      0       0       0
                        '4'      0       0       0
                        '5'      0       0       0
;

Table cHVL(index2,component)
                                 '2'     '1'     '3'
                        '1'      0       0       0
                        '2'      0       0       0
                        '3'      0       0       0
                        '4'      0       0       0
                        '5'      0       0       0
;

alias(component,component2)
set index3 /1*3/;
Table cGamma(index3,component,component2)
         '1'     '2'     '3'
'1'.'1'  0        0      0
'1'.'2'  0        0      0
'1'.'3'  0        0      0
'2'.'1'  0        0      0
'2'.'2'  0        0      0
'2'.'3'  0        0      0
'3'.'1'  0        0      0
'3'.'2'  0        0      0
'3'.'3'  0        0      0
;

$gdxout properties/AcetoneChloroformBenzene
$unload component, P, Tref, cPr, cHV, cHVL, cGamma
$gdxout
