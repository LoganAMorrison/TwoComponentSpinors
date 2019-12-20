(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

FermionSpinSum::usage="";
	
Begin["Private`"]


(* WeylLineApplyVectorFierz[expr]
	Reduce expression of the form:
		(z1.S_mu.z2d) * (z3d.Sb_mu.z4)
		(z1d.Sd_mu.z2) * (z3d.Sb_mu.z4)
		(z1.S_mu.z2d) * (z3.S_mu.z4d)
*)
WeylLineApplyVectorFierz[expr___]:=expr//.{
	(* (z1.sb_mu.z2) * (z3.s_mu.z4) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWMR[iLTS[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWML[iLTS[mu_]]] :>
		-2 * WeylLine[{s1,p1,m1},{s4,p4,m4},iWML[Weyl1]] *
		WeylLine[{s2,p2,m2},{s3,p3,m3},iWMR[Weyl1]],
	(* (z1.sb_mu.z2) * (z3.sb_mu.z4) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWML[iLTS[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWML[iLTS[mu_]]] :>
		2 * WeylLine[{s1,p1,m1},{s3,p3,m3},iWMR[Weyl1]] *
		WeylLine[{s4,p4,m4},{s2,p2,m2},iWML[Weyl1]],
	(* (z1.s_mu.z2) * (z3.s_mu.z4) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWML[iLTS[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWML[iLTS[mu_]]] :>
		2 * WeylLine[{s1,p1,m1},{s3,p3,m3},iWMR[Weyl1]] *
		WeylLine[{s4,p4,m4},{s2,p2,m2},iWML[Weyl1]]
};

(* WeylLinesRemoveZero[expr]
	Kill off factors of z1.z1 or z1d.z1d since z's are communting.
*)
WeylLinesRemoveZero[expr___]:=expr//.{
	WeylLine[{s1_,p1_,m1_}, {s1_,p1_,m1_}, (pat:iWML|iWMR)[Weyl1]] :> 0
};


(* WeylLineMetricContract[expr]
	Contracts adjacent WeylLines of the form:
		g_{mu,nu} * WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_nu]]
	into
		WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_mu]]
*)
WeylLineMetricContract[expr___]:=expr//.{
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[iLTS[mu1_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},(pat2:iWML|iWMR)[iLTS[mu2_]]] *
	LTensor[MetricG,mu1_,mu2_] :>
		WeylLine[{s1,p1,m1},{s2,p2,m2},pat1[iLTS[mu1]]] *
		WeylLine[{s3,p3,m3},{s4,p4,m4},pat2[iLTS[mu1]]]
};

(* WeylLineFixOrder[expr]
	Reverses the order of spinors which aren't facing tail to head, i.e.:
		(z_1.X.z_2) * (z_3.Y.z_2)
	get's converted into
		(z_1.X.z_2) * (z_2.Y'.z_3)
*)
WeylLineFixOrder[expr___]:=expr/.{
	(* (z1.1.z2)^2 -> -(z1.1.z2) * (z2.1.z1) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat:iWML|iWMR)[Weyl1]]^2 :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},pat[Weyl1]] *
		(* change iWML->iWMR and iWMR->iWML *)
		WeylLine[{s2,p2,m2},{s1,p1,m1},
			(Switch[pat,iWML,iWMR,iWMR,iWML])[Weyl1]
		],
	(* (z1.X.z2) (z3.1.z2) -> -(z1.X.z2) * (z2.1.z3) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[mtx___]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p2_,m2_},(pat2:iWML|iWMR)[Weyl1]] :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},pat1[mtx]] *
		WeylLine[{s4,p2,m2},{s3,p3,m3},
			(Switch[pat2,iWML,iWMR,iWMR,iWML])[Weyl1]
		],
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[mtx___]] *
	WeylLine[{s3_,p1_,m1_},{s4_,p4_,m4_},(pat2:iWML|iWMR)[Weyl1]] :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},pat1[mtx]] *
		WeylLine[{s4,p4,m4},{s3,p1,m1},
			(Switch[pat2,iWML,iWMR,iWMR,iWML])[Weyl1]
		],
	(* (z1.X.z2) (z3.Y.z2) -> (-1)^(Length[Y]+1) (z1.X.z2) * (z2.YR.z3)
		with YR being the reverse of Y with WeylS<->WeylSBar *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[mtx1___]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p2_,m2_},(pat2:iWML|iWMR)[mtx2___]] :>
		(*if even # of sigmas, we get extra minus sign *)
		(-1)^(Count[List[mtx2], WeylS]+1) *
		WeylLine[{s1,p1,m1},{s2,p2,m2},pat1[mtx1]] *
		(* if even # of sigmas, don't change head of pat2, else swap *)
		WeylLine[{s4,p2,m2},{s3,p3,m3},
			(If[EvenQ[Count[List[mtx2], WeylS]], pat2,
				Switch[pat2,iWML,iWMR,iWMR,iWML]
				])@@Reverse[List[mtx2]]
		],
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[mtx1___]] *
	WeylLine[{s3_,p1_,m1_},{s4_,p4_,m4_},(pat2:iWML|iWMR)[mtx2___]] :>
		(*if even # of sigmas, we get extra minus sign *)
		(-1)^(Count[List[mtx2], WeylS]+1) *
		WeylLine[{s1,p1,m1},{s2,p2,m2},pat1[mtx1]] *
		(* if even # of sigmas, don't change head of pat2, else swap *)
		WeylLine[{s4,p4,m4},{s3,p1,m1},
			(If[EvenQ[Count[List[mtx2], WeylS]], pat2,
				Switch[pat2,iWML,iWMR,iWMR,iWML]
				])@@Reverse[List[mtx2]]
		]
};

(* WeylLineCollapse[expr]
	Collapses two compatible adjacent WeylLines into one WeylLine. Compatible
	WeylLines are of the form:
		WeylLine[wfl,{s1,p,m},X] * WeylLine[{s2,p,m},wfr,Y]
	These are converted into
		WeylLine[wfl,wfr,X.Y]
*)
WeylLineCollapse[expr___]:= expr/.{
	(* X.xd.x.Y  or  X.x.xd.Y  or  X.y.yd.Y  or  X.yd.y.Y *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},(pat1:iWML|iWMR)[mtx1___]] *
	WeylLine[{s2_,p2_,m2_},{s4_,p4_,m4_},(pat2:iWML|iWMR)[mtx2___]] :>
		WeylLine[{s1,p1,m1},{s4,p4,m4},pat2[mtx1,LDot[p2, WeylS],mtx2]],

	(* X.x.y.Y  or  X.y.x.Y *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWML[mtx1___]] *
	WeylLine[{s3_,p2_,m2_},{s4_,p4_,m4_},(pat:iWML|iWMR)[mtx2___]] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},pat[mtx1,s2*m2*Weyl1,mtx2]
	]/;SameQ[s3,-s2],

	(* X.yd.xd.Y  or  X.xd.yd.Y *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWMR[mtx1___]] *
	WeylLine[{s3_,p2_,m2_},{s4_,p4_,m4_},(pat:iWML|iWMR)[mtx2___]] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},pat[mtx1,-s2*m2*Weyl1,mtx2]
	]/;SameQ[s3,-s2]
};


WeylLineToTrace[expr___]:=expr/.{
	(* x.x = y.y = xd.xd = yd.yd =0 *)
	WeylLine[{s_,p_,m_},{s_,p_,m_},(pat:iWML|iWMR)[Weyl1]]:> 0,

	(* xd.X.x  or  yd.X.y  or  x.X.xd  or  y.X.yd *)
	WeylLine[{s_,p_,m_},{s_,p_,m_}, (pat:iWML|iWMR)[mtx___]]:>
		WeylTrace[ConvertToExternal[
			(Switch[pat,iWML,iWMR,iWMR,iWML])[mtx,LDot[p,WeylS]]
		]],

	(* x.X.y  or  y.X.x *)
	WeylLine[{s1_,p_,m_},{s2_,p_,m_}, iWML[mtx___]]:>
		WeylTrace[ConvertToExternal[iWML[mtx,-s1*m*Weyl1]]]/;SameQ[s1,-s2],

	(* xd.X.yd  or  yd.X.xd *)
	WeylLine[{s1_,p_,m_},{s2_,p_,m_}, iWMR[mtx___]]:>
		WeylTrace[ConvertToExternal[iWMR[mtx,s1*m*Weyl1]]]/;SameQ[s1,-s2]
};


InternalFermionSpinSum[expr___]:=Module[{iexpr},
	iexpr=X`Utilities`Uncontract[#, WeylS]&/@expr;
	iexpr=WeylLineMetricContract[iexpr];
	iexpr=WeylLineApplyVectorFierz[iexpr];
	iexpr=WeylLinesRemoveZero[iexpr];
	iexpr=WeylLineFixOrder[iexpr];
	iexpr=FixedPoint[WeylLineCollapse,iexpr];
	iexpr=X`Utilities`Uncontract[#, WeylS]&/@iexpr;
	FixedPoint[WeylLineToTrace,iexpr]
];


FermionSpinSum[expr_]:=Module[{iexpr},
	iexpr=ExpandAll[ConvertToInternal[expr]];
	iexpr=FixedPoint[InternalFermionSpinSum,iexpr];
	iexpr(*ConvertToExternal[iexpr]*)
]

End[]
