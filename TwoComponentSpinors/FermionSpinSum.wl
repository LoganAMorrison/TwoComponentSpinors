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
	WeylLine[{s1_,p1_,m1_}, {s1_,p1_,m1_}, (pat:iWMLL|iWMRR)[Weyl1]] :> 0,
	WeylLine[{s1_,p1_,m1_}, {s1_,p1_,m1_}, (pat:iWMLL|iWMRR)[]] :> 0
};


(* WeylLineMetricContract[expr]
	Contracts adjacent WeylLines of the form:
		g_{mu,nu} * WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_nu]]
	into
		WeylLine[...,...,iWM[X_mu]] * WeylLine[...,...,iWM[X_mu]]
*)
WeylLineMetricContract[expr___]:=expr//.{
	WeylLine[kin1_List, kin2_List, (pat1:iWMLL|iWMLR|iWMRR|iWMRL)[left1___, iLTS[mu1_], right1___]] *
	WeylLine[kin3_List, kin4_List, (pat2:iWMLL|iWMLR|iWMRR|iWMRL)[left2___, iLTS[mu2_], right2___]] *
	LTensor[MetricG, mu1_, mu2_] :>
		WeylLine[kin1, kin2, pat1[left1, iLTS[mu1], right1]] *
		WeylLine[kin3, kin4, pat2[left2, iLTS[mu1], right2]]
};

(* WeylLineFixOrder[expr]
	Reverses the order of spinors which aren't facing tail to head, i.e.:
		(z_1.X.z_2) * (z_3.Y.z_2)
	get's converted into
		(z_1.X.z_2) * (z_2.Y'.z_3)
*)
WeylLineFixOrder[expr___]:=expr/.{
	(* (z1.X.z2) (z3.Y.z2) -> (-1)^(Length[Y]+1) (z1.X.z2) * (z2.YR.z3)
		with YR being the reverse of Y with WeylS<->WeylSBar *)
	WeylLineProduct[
		left___,
		WeylLine[{s1_, kin1__},{s2_, kin2__},
			(pat1:iWMLL|iWMRR|iWMLR|iWMRL)[mtx1___]],
		cent___,
		WeylLine[{s3_, kin3__},{s4_, kin2__},
			(pat2:iWMLL|iWMRR)[mtx2___]],
		right___
	] :> (-1) * (*if even # of sigmas, we get extra minus sign *)
		WeylLineProduct[
			left,cent,right,
			WeylLine[{s1,kin1},{s2,kin2}, pat1[mtx1]],
			WeylLine[{s4,kin2},{s3,kin3}, pat2@@Reverse[List[mtx2]]
			]
		],

	WeylLineProduct[
		left___,
		WeylLine[{s1_, kin1__}, {s2_, kin2__}, (pat1:iWMLL|iWMRR|iWMLR|iWMRL)[mtx1___]],
		cent___,
		WeylLine[{s3_, kin3__}, {s4_, kin2__}, (pat2:iWMLR|iWMRL)[mtx2___]],
		right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[{s1, kin1}, {s2, kin2}, pat1[mtx1]],
			(* swap iWMRL and iWMLR *)
			WeylLine[{s4, kin2}, {s3, kin3}, (Switch[pat2,iWMRL,iWMLR,iWMLR,iWMRL])@@Reverse[List[mtx2]]]
		(* Only flip if necessary. *)
		]
};

(* WeylLineCollapse[expr]
	Collapses two compatible adjacent WeylLines into one WeylLine. Compatible
	WeylLines are of the form:
		WeylLine[wfl,{s1,p,m},X] * WeylLine[{s2,p,m},wfr,Y]
	These are converted into
		WeylLine[wfl,wfr,X.Y]
*)
WeylLineCollapse[expr___]:= expr//.{
	(* (X.xd)*(x.Y)  or  (X.yd)*(y.Y) *)
	WeylLineProduct[
		left___, WeylLine[kin1_List, {s_,p_,m_}, iWMRR[mtx1___]],
		cent___, WeylLine[{s_,p_,m_}, kin3_List, iWMLL[mtx2___]],
		right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1, kin3, iWMRL[mtx1, LDot[p, WeylS], mtx2]]
	],

	WeylLineProduct[
		left___, WeylLine[kin1_List, {s_,p_,m_}, iWMRR[mtx1___]],
		cent___, WeylLine[{s_,p_,m_}, kin3_List, iWMLR[mtx2___]],
		right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1, kin3, iWMRR[mtx1, LDot[p, WeylS],mtx2]]
		],

	WeylLineProduct[
		left___, WeylLine[kin1_List, {s_,p_,m_}, iWMLR[mtx1___]],
		cent___, WeylLine[{s_,p_,m_}, kin3_List, iWMLL[mtx2___]],
		right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1, kin3, iWMLL[mtx1, LDot[p, WeylS], mtx2]]
		],

	WeylLineProduct[
		left___, WeylLine[kin1_List,{s_,p_,m_}, iWMLR[mtx1___]],
		cent___, WeylLine[{s_,p_,m_},kin3_List, iWMLR[mtx2___]],
		right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1, kin3, iWMLR[mtx1, LDot[p, WeylS], mtx2]]
		],

	(* X.x.xd.Y  or  X.y.yd.Y *)
	WeylLineProduct[
		left___, WeylLine[kin1_List, {s_,p_,m_}, iWMLL[mtx1___]],
		cent___, WeylLine[{s_,p_,m_}, kin3_List, iWMRL[mtx2___]], right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1, kin3, iWMLL[mtx1,LDot[p,WeylS],mtx2]]
		],

	WeylLineProduct[
		left___, WeylLine[kin1_List,{s_,p_,m_},iWMLL[mtx1___]],
		cent___, WeylLine[{s_,p_,m_},kin3_List,iWMRR[mtx2___]], right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1,kin3, iWMLR[mtx1,LDot[p,WeylS],mtx2]]
		],

	WeylLineProduct[
		left___, WeylLine[kin1_List,{s_,p_,m_},iWMRL[mtx1___]],
		cent___, WeylLine[{s_,p_,m_},kin3_List,iWMRL[mtx2___]], right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1,kin3, iWMRL[mtx1,LDot[p,WeylS],mtx2]]
		],

	WeylLineProduct[
		left___, WeylLine[kin1_List,{s_,p_,m_},iWMRL[mtx1___]],
		cent___, WeylLine[{s_,p_,m_},kin3_List,iWMRR[mtx2___]], right___
	] :> WeylLineProduct[left,cent,right,
			WeylLine[kin1,kin3, iWMRR[mtx1,LDot[p,WeylS],mtx2]]
		],


	(* X.x.y.Y  or  X.y.x.Y *)
	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMLL[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMLL[mtx2___]],
		right___
	] :> s2*m2*WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMLL[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMLL[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMLR[mtx2___]],
		right___
	] :> s2*m2*WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMLR[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMRL[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMLL[mtx2___]],
		right___
	] :> s2*m2*WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMRL[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List, {s2_,p2_,m2_}, iWMRL[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_}, kin4_List, iWMLR[mtx2___]],
		right___
	] :> s2*m2*WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMRR[mtx1,mtx2]]
	]/;SameQ[s3,-s2],


	(* XR.yd.xd.YL  or  XR.xd.yd.YL *)
	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMLR[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMRL[mtx2___]],
		right___
	] :> -s2*m2*WeylLineProduct[
		left,cent,right,
		WeylLine[kin1,kin4,iWMLL[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMLR[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMRR[mtx2___]],
		right___
	] :> -s2*m2*WeylLineProduct[
		left,cent,right,
		WeylLine[kin1,kin4,iWMLR[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMRR[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMRL[mtx2___]],
		right___
	] :> -s2*m2* WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMRL[mtx1,mtx2]]
	]/;SameQ[s3,-s2],

	WeylLineProduct[
		left___,WeylLine[kin1_List,{s2_,p2_,m2_},iWMRR[mtx1___]],
		cent___,WeylLine[{s3_,p2_,m2_},kin4_List,iWMRR[mtx2___]],
		right___
	] :> -s2*m2*WeylLineProduct[left,cent,right,
		WeylLine[kin1,kin4,iWMRR[mtx1,mtx2]]
	]/;SameQ[s3,-s2]
};


(* WeylLineToTrace[expr]
	Converts we WeylLine with same particle at endpoint into a
	trace. i.e.,
		WeylLine[{s1_,p_,m_},{s2_,p_,m_},X]
	gets converted into the a trace over X times the outer product
	of the external wavefunctions.
*)
WeylLineToTrace[expr___]:=expr//.{
	(* x.x = y.y = xd.xd = yd.yd =0 *)
	WeylLineProduct[WeylLine[{s_,p_,m_},{s_,p_,m_},(pat:iWMLL|iWMRR)[]]]:> 0,

	(* xd.X.x  or  yd.X.y  or  x.X.xd  or  y.X.yd *)
	WeylLineProduct[
		left___,
		WeylLine[{s_,p_,m_},{s_,p_,m_}, (pat:iWMRL|iWMLR)[mtx___]],
		right___
	]:> WeylLineProduct[left,right] *
			WeylTrace[
			(Switch[pat,iWMRL,iWMLL,iWMLR,iWMRR])[LDot[p,WeylS],mtx]
	],

	(* x.X.y  or  y.X.x or xd.X.yd  or  yd.X.xd *)
	WeylLineProduct[
		left___,
		WeylLine[{s1_,p_,m_},{s2_,p_,m_}, (pat:iWMLL|iWMRR)[mtx___]],
		right___
	]:> WeylLineProduct[left,right] * Switch[pat,iWMLL,s2,iWMRR,s1] * m *
			WeylTrace[pat[mtx]
	]/;SameQ[s1,-s2]
};


InternalFermionSpinSum[WeylLineProduct[lines___]]:=Module[{iexpr},
	iexpr=WeylLineProduct[lines];
	(*
	iexpr=X`Utilities`Uncontract[#, WeylS]&/@iexpr;
	iexpr=WeylLineMetricContract[iexpr];
	iexpr=WeylLineApplyVectorFierz[iexpr];
	iexpr=WeylLinesRemoveZero[iexpr];
	*)
	(*iexpr=WeylLineFixOrder[iexpr];*)
	iexpr=FixedPoint[WeylLineFixOrder[WeylLineCollapse[#]]&,iexpr];
	Contract[WeylLineToTrace[iexpr]]
];


FermionSpinSum[expr_]:=Module[{iexpr},
	iexpr=ConvertToInternal[Collect[
		Expand[expr]/.{
			wlp_WeylLineProduct:>
				WeylLineProductExpand[wlp]
		},
		WeylLineProduct[___]]
	];
	iexpr=iexpr/.{lines_WeylLineProduct:>InternalFermionSpinSum[lines]};
	ConvertToExternal[iexpr]
]

End[]
