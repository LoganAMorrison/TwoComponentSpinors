(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

FermionSpinSum::usage="";
WeylLineCollapse::usage="";
ConvertToInternal::usage="";
WeylLineFixOrder::usage="";

Begin["Private`"]


(* ConvertToInternal[expr]
	Convert expressions into more convinient form.
*)
ConvertToInternal[expr___]:=expr//.{
	WeylMatrix[mtx___]:>iWM[mtx],
	LTensor[WeylS, mu_] :> iLTS[mu],
	LTensor[WeylSBar, mu_] :> iLTSB[mu]
};

(* ConvertToInternal[expr]
	Convert expressions back to external form.
*)
ConvertToExternal[expr___]:=expr//.{
	iWM[mtx___] :> WeylMatrix[mtx],
	iLTS[mu_] :> LTensor[WeylS, mu],
	iLTSB[mu_] :> LTensor[WeylSBar, mu]
};


(* WeylLineApplyVectorFierz[expr]
	Reduce expression of the form:
		(z1.S_mu.z2d) * (z3d.Sb_mu.z4)
		(z1d.Sd_mu.z2) * (z3d.Sb_mu.z4)
		(z1.S_mu.z2d) * (z3.S_mu.z4d)
*)
WeylLineApplyVectorFierz[expr___]:=expr//.{
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[iLTS[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWM[iLTSB[mu_]]] :>
		-2 * WeylLine[{s1,p1,m1},{s4,p4,m4},iWM[Weyl1]] *
		WeylLine[{s2,p2,m2},{s3,p3,m3},iWM[Weyl1Bar]],

	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[iLTSB[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWM[iLTSB[mu_]]] :>
		2 * WeylLine[{s1,p1,m1},{s3,p3,m3},iWM[Weyl1Bar]] *
		WeylLine[{s4,p4,m4},{s2,p2,m2},iWM[Weyl1]],

	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[iLTS[mu_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWM[iLTS[mu_]]] :>
		2 * WeylLine[{s1,p1,m1},{s3,p3,m3},iWM[Weyl1]] *
		WeylLine[{s4,p4,m4},{s2,p2,m2},iWM[Weyl1Bar]]
};

(* WeylLinesRemoveZero[expr]
	Kill off factors of z1.z1 or z1d.z1d since z's are communting.
*)
WeylLinesRemoveZero[expr___]:=expr//.{
	WeylLine[{s1_,p1_,m1_}, {s1_,p1_,m1_}, iWM[Weyl1|Weyl1Bar]] :> 0
};


(* WeylLineMetricContract[expr]
	Contracts adjacent WeylLines of the form:
		g_{mu,nu} * WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_nu]]
	into
		WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_mu]]
*)
WeylLineMetricContract[expr___]:=expr//.{
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[(pat1:iLTS|iLTSB)[mu1_]]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p4_,m4_},iWM[(pat2:iLTS|iLTSB)[mu2_]]] *
	LTensor[MetricG,mu1_,mu2_] :>
		WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[pat1[mu1]]] *
		WeylLine[{s3,p3,m3},{s4,p4,m4},iWM[pat2[mu1]]]
};

(* WeylLineFixOrder[expr]
	Reverses the order of spinors which aren't facing tail to head, i.e.:
		(z_1.X.z_2) * (z_3.Y.z_2)
	get's converted into
		(z_1.X.z_2) * (z_2.Y'.z_3)
*)
WeylLineFixOrder[expr___]:=ReleaseHold[expr/.{
	(* (z1.1.z2)^2 -> -(z1.1.z2) * (z2.1.z1) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[pat1:Weyl1|Weyl1Bar]]^2 :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[pat1]] *
		WeylLine[{s2,p2,m2},{s1,p1,m1},iWM[pat1]],
	(* (z1.X.z2) (z3.1.z2) -> -(z1.X.z2) * (z2.1.z3) *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[mtx___]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p2_,m2_},iWM[pat2:Weyl1|Weyl1Bar]] :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[mtx]] *
		WeylLine[{s4,p2,m2},{s3,p3,m3},iWM[pat2]],
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[mtx___]] *
	WeylLine[{s3_,p1_,m1_},{s4_,p4_,m4_},iWM[pat2:Weyl1|Weyl1Bar]] :>
		-1 * WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[mtx]] *
		WeylLine[{s4,p4,m4},{s3,p1,m1},iWM[pat2]],
	(* (z1.X.z2) (z3.Y.z2) -> (-1)^(Length[Y]+1) (z1.X.z2) * (z2.YR.z3)
		with YR being the reverse of Y with WeylS<->WeylSBar *)
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[mtx1___]] *
	WeylLine[{s3_,p3_,m3_},{s4_,p2_,m2_},iWM[mtx2___]] :>
		(-1)^(Length[List[mtx2]]+1) *
		WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[mtx1]] *
			Hold[ReplaceAll][
				WeylLine[{s4,p2,m2},{s3,p3,m3},iWM@@Reverse[List[mtx2]]
				], {iLTS->iLTSB,iLTSB->iLTS}],
	WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},iWM[mtx1___]] *
	WeylLine[{s3_,p1_,m1_},{s4_,p4_,m4_},iWM[mtx2___]] :>
		(-1)^(Length[List[mtx2]]+1) *
		WeylLine[{s1,p1,m1},{s2,p2,m2},iWM[mtx1]] *
			Hold[ReplaceAll][
				WeylLine[{s4,p4,m4},{s3,p1,m1},iWM@@Reverse[List[mtx2]]
				], {iLTS->iLTSB,iLTSB->iLTS}]
}]

(* WeylLineCollapse[expr]
	Collapses two compatible adjacent WeylLines into one WeylLine. Compatible
	WeylLines are of the form:
		WeylLine[wfl,{s1,p,m},X] * WeylLine[{s2,p,m},wfr,Y]
	These are converted into
		WeylLine[wfl,wfr,X.Y]
*)
WeylLineCollapse[expr___]:= expr/.{
	(* ...X xd x Y... *)
	WeylLine[{s1_,p1_,m1_},{1,p2_,m2_},
		iWM[left___,pat1:Weyl1Bar|iLTS[mu_]]
	] *
	WeylLine[{1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1|iLTS[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,LDot[p2, WeylSBar],pat2,right]
	],
	(* ...X x xd Y... *)
	WeylLine[{s1_,p1_,m1_},{1,p2_,m2_},
		iWM[left___,pat1:Weyl1|iLTSB[mu_]]
	] *
	WeylLine[{1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1Bar|iLTSB[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,LDot[p2, WeylS],pat2,right]
	],


	(* ...X yd y Y... *)
	WeylLine[{s1_,p1_,m1_},{-1,p2_,m2_},
		iWM[left___,pat1:Weyl1Bar|iLTS[mu_]]
	] *
	WeylLine[{-1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1|iLTS[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,LDot[p2, WeylSBar],pat2,right]
	],
	(* ...X y yd Y... *)
	WeylLine[{s1_,p1_,m1_},{-1,p2_,m2_},
		iWM[left___,pat1:Weyl1|iLTSB[mu_]]
	] *
	WeylLine[{-1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1Bar|iLTSB[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,LDot[p2, WeylS],pat2,right]
	],


	(* ...X x y Y... *)
	WeylLine[{s1_,p1_,m1_},{1,p2_,m2_},
		iWM[left___,pat1:Weyl1|iLTSB[mu_]]
	] *
	WeylLine[{-1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1|iLTS[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,m2*Weyl1,pat2,right]
	],
	(* ...X y x Y... *)
	WeylLine[{s1_,p1_,m1_},{-1,p2_,m2_},
		iWM[left___,pat1:Weyl1|iLTSB[mu_]]
	] *
	WeylLine[{1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1|iLTS[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,-m2*Weyl1,pat2,right]
	],


	(* ...X yd xd Y... *)
	WeylLine[{s1_,p1_,m1_},{-1,p2_,m2_},
		iWM[left___,pat1:Weyl1Bar|iLTS[mu_]]
	] *
	WeylLine[{1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1Bar|iLTSB[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,m2*Weyl1Bar,pat2,right]
	],
	(* ...X xd yd Y... *)
	WeylLine[{s1_,p1_,m1_},{1,p2_,m2_},
		iWM[left___,pat1:Weyl1Bar|iLTS[mu_]]
	] *
	WeylLine[{-1,p2_,m2_},{s4_,p4_,m4_},
		iWM[pat2:Weyl1Bar|iLTSB[nu_],right___]
	] :>
	WeylLine[{s1,p1,m1},{s4,p4,m4},
		iWM[left,pat1,-m2*Weyl1Bar,pat2,right]
	]
};


WeylLineToTrace[expr___]:=expr/.{
	(* x.x =0 *)
	WeylLine[{1,p_,m_},{1,p_,m_},iWM[Weyl1]]:> 0,
	(* y.y =0 *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},iWM[Weyl1]]:> 0,
	(* x.y *)
	WeylLine[{1,p_,m_},{-1,p_,m_},iWM[Weyl1]]:> 2*m,
	(* y.x *)
	WeylLine[{-1,p_,m_},{1,p_,m_},iWM[Weyl1]]:> -2*m,

	(* xd.xd =0 *)
	WeylLine[{1,p_,m_},{1,p_,m_},iWM[Weyl1Bar]]:> 0,
	(* yd.yd =0 *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},iWM[Weyl1Bar]]:> 0,
	(* xd.yd *)
	WeylLine[{1,p_,m_},{-1,p_,m_},iWM[Weyl1Bar]]:> -2*m,
	(* yd.xd *)
	WeylLine[{-1,p_,m_},{1,p_,m_},iWM[Weyl1Bar]]:> 2*m,

	(* xd.Sb_mu.x *)
	WeylLine[{1,p_,m_},{1,p_,m_},iWM[iLTSB[mu_]]]:>
		WeylTrace[ConvertToExternal[iWM[iLTSB[mu],LDot[p,WeylS]]]],
	(* x.S_mu.xd *)
	WeylLine[{1,p_,m_},{1,p_,m_},iWM[iLTS[mu_]]]:>
		WeylTrace[ConvertToExternal[iWM[iLTS[mu],LDot[p,WeylSBar]]]],
	(* yd.Sb_mu.y *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},iWM[iLTSB[mu_]]]:>
		WeylTrace[ConvertToExternal[iWM[iLTSB[mu],LDot[p,WeylS]]]],
	(* y.S_mu.yd *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},iWM[iLTS[mu_]]]:>
		WeylTrace[ConvertToExternal[iWM[iLTS[mu],LDot[p,WeylSBar]]]],


	(* x.X.y *)
	WeylLine[{1,p_,m_},{-1,p_,m_},
		iWM[pat1:iLTS[mu_]|Weyl1,cent__,pat2:iLTSB[nu_]|Weyl1]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,-m*Weyl1]]],
	(* xd.X.yd *)
	WeylLine[{1,p_,m_},{-1,p_,m_},
		iWM[pat1:iLTSB[mu_]|Weyl1Bar,cent__,pat2:iLTS[nu_]|Weyl1Bar]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,m*Weyl1Bar]]],
	(* y.X.x *)
	WeylLine[{-1,p_,m_},{1,p_,m_},
		iWM[pat1:iLTS[mu_]|Weyl1,cent__,pat2:iLTSB[nu_]|Weyl1]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,m*Weyl1]]],
	(* yd.X.xd *)
	WeylLine[{-1,p_,m_},{1,p_,m_},
		iWM[pat1:iLTSB[mu_]|Weyl1Bar,cent__,pat2:iLTS[nu_]|Weyl1Bar]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,-m*Weyl1Bar]]],

	(* xd.X.x *)
	WeylLine[{1,p_,m_},{1,p_,m_},
		iWM[pat1:iLTSB[mu_]|Weyl1Bar,cent__,pat2:iLTSB[nu_]|Weyl1]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,LDot[p,WeylS]]]],
	(* yd.X.y *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},
		iWM[pat1:iLTSB[mu_]|Weyl1Bar,cent__,pat2:iLTSB[nu_]|Weyl1]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,LDot[p,WeylS]]]],
	(* x.X.xd *)
	WeylLine[{1,p_,m_},{1,p_,m_},
		iWM[pat1:iLTS[mu_]|Weyl1,cent__,pat2:iLTS[nu_]|Weyl1Bar]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,LDot[p,WeylSBar]]]],
	(* y.X.yd *)
	WeylLine[{-1,p_,m_},{-1,p_,m_},
		iWM[pat1:iLTS[mu_]|Weyl1,cent__,pat2:iLTS[nu_]|Weyl1Bar]]:>
		WeylTrace[ConvertToExternal[iWM[pat1,cent,pat2,LDot[p,WeylSBar]]]]
};


InternalFermionSpinSum[expr___]:=Module[{iexpr},
	iexpr=X`Utilities`Uncontract[#, WeylS|WeylSBar]&/@expr;
	iexpr=WeylLineMetricContract[iexpr];
	iexpr=WeylLineApplyVectorFierz[iexpr];
	iexpr=WeylLinesRemoveZero[iexpr];
	iexpr=WeylLineFixOrder[iexpr];
	Print[];
	Print[];
	Print[iexpr];
	Print[];
	Print[];
	iexpr=FixedPoint[WeylLineCollapse,iexpr];
	iexpr=X`Utilities`Uncontract[#, WeylS|WeylSBar]&/@iexpr;
	(*Print[iexpr];*)
	FixedPoint[WeylLineToTrace,iexpr]
];


FermionSpinSum[expr_]:=Module[{iexpr},
	(*
	iexpr=PreprocessWeylProducts[ExpandAll[expr]];

	iexpr/.{
		WeylLine[{s1_,p1_,m1_},{s2_,p2_,m2_},WeylMatrix[mtx1___]] *
		WeylLine[{s3_,p2_,m2_},{s4_,p1_,m1_}, WeylMatrix[mtx2___]] :>
		InternalFermionSpinSum[
			WeylLine[{s1,p1,m1},{s2,p2,m2},WeylMatrix[mtx1]],
			WeylLine[{s3,p2,m2},{s4,p1,m1},WeylMatrix[mtx2]]
		]
	}*)
	iexpr=ExpandAll[ConvertToInternal[expr]];
	iexpr=FixedPoint[InternalFermionSpinSum,iexpr];
	ConvertToExternal[iexpr]
]

End[]
