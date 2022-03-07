(* ::Package:: *)

(* ::Section::Initialization:: *)
(*Initialization Group*)


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]]


(* ::Input::Initialization:: *)
$Assumptions = {c>0}


(* ::Input::Initialization:: *)
RulePN = {\[Mu]-> mu, \[Nu]-> nu}


(* ::Input::Initialization:: *)
RulePM = {E1[p]-> E1, E2[p]-> E2, EE[p]-> EE, \[Gamma][p]-> gamma, \[Sigma][p]-> sig,\[Xi][p]-> xi, \[Nu]-> nu}


(* ::Input::Initialization:: *)
RuleNotation = {Dot[1, y_]->y,Dot[y_,1]-> y}


(* ::Input::Initialization:: *)
RuleDerivs = {Derivative[1][E1][p]-> dE1, Derivative[1][E2][p]-> dE2,Derivative[1][\[Gamma]][p]-> 1/(m c^2)(dE1+dE2), Derivative[1][\[Xi]][p]-> 1/(EE^3)(E1-E2)(E1 dE2 - E2 dE1), Derivative[1][\[Sigma]][p]-> 1/(m1 m2 c^2)(E1/c^2 dE2 + E2/c^2 dE1 + 2 p)}


(* ::Input::Initialization:: *)
pr = (p . r)/Sqrt[(r . r)]


(* ::Input::Initialization:: *)
H0PN = (p . p)/2 - 1/Sqrt[(r . r)]


(* ::Input::Initialization:: *)
H0PNfull = m c^2 + \[Mu] H0PN


(* ::Input::Initialization:: *)
E0PN = FortranForm[H0PNfull/.RulePN/.RuleNotation]


(* ::Input::Initialization:: *)
Export["E0PN.txt",{E0PN}]


(* ::Input::Initialization:: *)
dHdr0PN = FortranForm[D[H0PN,r]/.RulePN/.RuleNotation]


(* ::Input::Initialization:: *)
Export["dHdr0PN.txt",{dHdr0PN}]


(* ::Input::Initialization:: *)
dHdp0PN =FortranForm[D[H0PN,p]/.RulePN/.RuleNotation]


(* ::Input::Initialization:: *)
Export["dHdp0PN.txt",{dHdp0PN}]


(* ::Input::Initialization:: *)
H1PN = H0PN + (1/8(3 \[Nu] -1)(p . p)^2 -1/2((3+\[Nu])(p . p)+\[Nu] pr^2)/Sqrt[(r . r)]+1/(2 (r . r)))/c^2


(* ::Input::Initialization:: *)
H1PNfull = m c^2 + \[Mu] H1PN


(* ::Input::Initialization:: *)
E1PN = FortranForm[H1PNfull/.RulePN/.RuleNotation]


(* ::Input::Initialization:: *)
Export["E1PN.txt",{E1PN}]


(* ::Input::Initialization:: *)
D[H1PN,r]/.RuleNotation


(* ::Input::Initialization:: *)
dHdr1PN = FortranForm[D[H1PN,r]/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdr1PN.txt",{dHdr1PN}]


(* ::Input::Initialization:: *)
D[H1PN,p]/.RuleNotation


(* ::Input::Initialization:: *)
dHdp1PN =FortranForm[D[H1PN,p]/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdp1PN.txt",{dHdp1PN}]


(* ::Input::Initialization:: *)
H2PN = H1PN + (1/16(1-5 \[Nu] + 5 \[Nu]^2)(p . p)^3 + 1/8((5 -20 \[Nu] - 3 \[Nu]^2)(p . p)^2 - 2 \[Nu]^2 pr^2 (p . p) - 3 \[Nu]^2 pr^4)/Sqrt[(r . r)]+1/2((5+ 8 \[Nu])(p . p) + 3 \[Nu] pr^2)/(r . r) - 1/4(1+3 \[Nu])/Sqrt[(r . r)]^3)/c^4//Simplify


(* ::Input::Initialization:: *)
H2PNfull = m c^2 + \[Mu] H2PN


(* ::Input::Initialization:: *)
E2PN = FortranForm[H2PNfull/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["E2PN.txt",{E2PN}]


(* ::Input::Initialization:: *)
dHdr2PN = FortranForm[D[H2PN,r]/.RuleNotation/.RulePN//Simplify]


(* ::Input::Initialization:: *)
Export["dHdr2PN.txt",{dHdr2PN}]


(* ::Input::Initialization:: *)
dHdp2PN =FortranForm[D[H2PN,p]/.RuleNotation/.RulePN//Simplify]


(* ::Input::Initialization:: *)
Export["dHdp2PN.txt",{dHdp2PN}]


(* ::Input::Initialization:: *)
H3PN = H2PN + (1/128(-5+35 \[Nu] -70 \[Nu]^2 + 35 \[Nu]^3)(p . p)^4 +1/16((-7+42 \[Nu] -53 \[Nu]^2 - 5 \[Nu]^3)(p . p)^3 + (2-3 \[Nu])\[Nu]^2 pr^2 (p . p)^2 + 3(1-\[Nu])\[Nu]^2 pr^4(p . p) - 5 \[Nu]^3 pr^6)/Sqrt[(r . r)]+(1/16(-27 + 136 \[Nu] + 109 \[Nu]^2)(p . p)^2 + 1/16(17 + 30 \[Nu])\[Nu] pr^2 (p . p) + 1/12(5 + 43 \[Nu])\[Nu] pr^4)/(r . r) +((-25/8 + (1/64 pi^2 - 335/48)\[Nu] - 23/8 \[Nu]^2)(p . p)+(-85/16-3/64 pi^2 -7/4 \[Nu])\[Nu] pr^2)/Sqrt[(r . r)]^3 + (1/8+(109/1
2 - 21/32 pi^2)\[Nu])/(r . r)^2)/c^6//Simplify


(* ::Input::Initialization:: *)
H3PNfull = m c^2 + \[Mu] H3PN


(* ::Input::Initialization:: *)
E3PN = FortranForm[H3PNfull/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["E3PN.txt",{E3PN}]


(* ::Input::Initialization:: *)
dHdr3PN = FortranForm[D[H3PN,r]/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdr3PN.txt",{dHdr3PN}]


(* ::Input::Initialization:: *)
dHdp3PN =FortranForm[D[H3PN,p]/.RulePN/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdp3PN.txt",{dHdp3PN}]


(* ::Input::Initialization:: *)
RulePM = {E1[p]-> E1, E2[p]-> E2, EE[p]-> EE, \[Gamma][p]-> gamma, \[Sigma][p]-> sig,\[Xi][p]-> xi, \[Nu]-> nu}


(* ::Input::Initialization:: *)
c1 = \[Nu]^2 m^2/(\[Gamma][p]^2 \[Xi][p])(1- 2 \[Sigma][p]^2)


(* ::Input::Initialization:: *)
H1PM =E1[p] + E2[p] + G/Sqrt[(r . r)] c1


(* ::Input::Initialization:: *)
E1PM = FortranForm[H1PM/.RulePM/.RuleNotation]


(* ::Input::Initialization:: *)
Export["E1PM.txt",{E1PM}]


(* ::Input::Initialization:: *)
dHdr1PM = FortranForm[D[H1PM,r]/.RulePM/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdr1PM.txt",{dHdr1PM}]


(* ::Input::Initialization:: *)
RuleDerivs = {Derivative[1][E1][p]-> dE1, Derivative[1][E2][p]-> dE2,Derivative[1][\[Gamma]][p]-> 1/(m c^2)(dE1+dE2), Derivative[1][\[Xi]][p]-> 1/(EE^3)(E1-E2)(E1 dE2 - E2 dE1), Derivative[1][\[Sigma]][p]-> 1/(m1 m2 c^2)(E1/c^2 dE2 + E2/c^2 dE1 + 2 p)}


(* ::Input::Initialization:: *)
dHdp1PM =FortranForm[D[H1PM,p]/.RuleNotation/.RuleDerivs/.RulePM//Simplify]


(* ::Input::Initialization:: *)
Export["dHdp1PM.txt",{dHdp1PM}]


(* ::Input::Initialization:: *)
c2 = \[Nu]^2 m^3/(\[Gamma][p]^2 \[Xi][p])(3/4(1-5 \[Sigma][p]^2)-4 \[Nu] \[Sigma][p](1- 2 \[Sigma][p]^2)/(\[Gamma][p] \[Xi][p])-\[Nu]^2(1-\[Xi][p])(1- 2 \[Sigma][p]^2)^2/(2 \[Gamma][p]^3 \[Xi][p]^2))


(* ::Input::Initialization:: *)
H2PM =E1[p] + E2[p] + G/Sqrt[(r . r)] c1+G^2/(r . r)/c^2c2


(* ::Input::Initialization:: *)
E2PM = FortranForm[H2PM/.RulePM/.RuleNotation]


(* ::Input::Initialization:: *)
Export["E2PM.txt",{E2PM}]


(* ::Input::Initialization:: *)
dHdr2PM = FortranForm[D[H2PM,r]/.RulePM/.RuleNotation//Simplify]


(* ::Input::Initialization:: *)
Export["dHdr2PM.txt",{dHdr2PM}]


(* ::Input::Initialization:: *)
dHdp2PM =FortranForm[D[H2PM,p]/.RuleNotation/.RuleDerivs/.RulePM//Simplify]


(* ::Input::Initialization:: *)
Export["dHdp2PM.txt",{dHdp2PM}]


(* ::Input::Initialization:: *)
c3 = \[Nu]^2 m^4/(\[Gamma][p]^2 \[Xi][p])(1/12(3-6 \[Nu] + 206 \[Nu] \[Sigma][p]-54 \[Sigma][p]^2+108 \[Nu] \[Sigma][p]^2 + 4 \[Nu] \[Sigma][p]^3)-4 \[Nu](3 + 12 \[Sigma][p]^2 - 4 \[Sigma][p]^4) ArcSinh[(\[Sigma][p]-1)/2]/Sqrt[\[Sigma][p]^2 -1]-3 \[Nu] \[Gamma][p](1-2 \[Sigma][p]^2)(1-5 \[Sigma][p]^2)/(2(1+\[Gamma][p])(1+\[Sigma][p]))-3 \[Nu] \[Sigma][p](7- 20 \[Sigma][p]^2)/(2 \[Gamma][p] \[Xi][p])+2 \[Nu]^3(3-4 \[Xi][p])\[Sigma][p](1-2 \[Sigma][p]^2)^2/(\[Gamma][p]^4 \[Xi][p]^3)-\[Nu]^2(3+8 \[Gamma][p]-3 \[Xi][p]-15 \[Sigma][p]^2-80 \[Gamma][p]\[Sigma][p]^2+15 \[Xi][p]\[Sigma][p]^2)(1-2 \[Sigma][p]^2)/(4 \[Gamma][p]^3 \[Xi][p]^2)+ \[Nu]^4(1-2 \[Xi][p])(1-2 \[Sigma][p]^2)^3/(2 \[Gamma][p]^6 \[Xi][p]^4))


(* ::Input::Initialization:: *)
H3PMonly = G^3/(Sqrt[r . r]^3)/c^4 c3


(* ::Input::Initialization:: *)
H3PM = H2PM + H3PMonly


(* ::Input::Initialization:: *)
E3PM = FortranForm[H3PM/.RulePM/.RuleNotation]


(* ::Input::Initialization:: *)
Export["E3PM.txt",{E3PM}]


(* ::Input::Initialization:: *)
dHdr3PMtemp = FortranForm[D[H3PMonly,r]/.RuleNotation/.RulePM//Simplify]


(* ::Input::Initialization:: *)
dHdr3PM =dHdr3PMtemp + dHdr2PM 


(* ::Input::Initialization:: *)
Export["dHdr3PM.txt",{dHdr3PM}]


(* ::Input::Initialization:: *)
dHdp3PMtemp = FortranForm[D[H3PMonly,p]/.RuleNotation/.RuleDerivs/.RulePM//Simplify]


(* ::Input::Initialization:: *)
dHdp3PM = dHdp3PMtemp + dHdp2PM 


(* ::Input::Initialization:: *)
Export["dHdp3PM.txt",{dHdp3PM}]
