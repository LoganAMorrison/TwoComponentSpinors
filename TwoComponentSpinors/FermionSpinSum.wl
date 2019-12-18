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


ConvertWeylLine[
	WeylLine[
		wf1_List /; SameQ[Length[wf1], 3],
		wf2_List /; SameQ[Length[wf2], 3], WeylMatrix[mtx___]
	]
] := Module[{wl,WM,LTS,LTSB,WMSS,WMSSB,WMSBS,WMSBSB},
	wl=WeylLine[wf1, wf2, WM[mtx]];

	wl=wl/.{LTensor[WeylS, mu_]:>LTS[mu]};
	wl=wl/.{LTensor[WeylSBar, mu_]:>LTSB[mu]};
	wl=wl/.{WeylMatrix[a___]:>WM[a]};

	(* Eliminate any inner WeylMatrix *)
	wl=wl//.{WM[a___,WM[b___],c___]:>WM[a,b,c]};

	wl=wl/.{
		WM[c1_. Weyl1 + a1_.] :> WMSS[c1 * Weyl1 + a1],
		WM[c1_. Weyl1Bar + a1_.] :> WMSBSB[c1 * Weyl1Bar + a1]
	};

	wl=wl//.{
		(* Replace WM[WeylS] with WMSS[WeylS] *)
		WM[c1_. LTS[mu1_] + a1_.] :> WMSS[c1 * LTS[mu1] + a1],
		(* Replace WM[WeylSBar] with WMSBSB[WeylSBar] *)
		WM[c1_. LTSB[mu1_] + a1_.] :> WMSBSB[c1 * LTSB[mu1] + a1],
		(* Replace WM[WeylS,...,WeylS] with WMSS[WeylS,...,WeylS] *)
		WM[c1_. LTS[mu1_] + a1_., a2___, c2_. LTS[mu2_] + a3_.] :>
			WMSS[c1 * LTS[mu1] + a1, a2, c2 * LTS[mu2] + a3],
		(* Replace WM[WeylS,...,WeylSBar] with WMSSB[WeylS,...,WeylSBar] *)
		WM[c1_. LTS[mu1_] + a1_., a2___, c2_. LTSB[mu2_] + a3_.] :>
			WMSSB[c1 * LTS[mu1] + a1, a2, c2 * LTSB[mu2] + a3],
		(* Replace WM[WeylSBar,...,WeylS] with WMSBS[WeylSBar,...,WeylS] *)
		WM[c1_. LTSB[mu1_] + a1_., a2___, c2_. LTS[mu2_] + a3_.] :>
			WMSBS[c1 * LTSB[mu1] + a1, a2, c2 * LTS[mu2] + a3],
		(* Replace WM[WeylSBar,...,WeylSBar] with WMSBSB[WeylSB,...,WeylSBar] *)
		WM[c1_. LTSB[mu1_] + a1_., a2___, c2_. LTSB[mu2_] + a3_.] :>
			WMSBSB[c1 * LTSB[mu1] + a1, a2, c2 * LTSB[mu2] + a3]
	};

	wl=wl/.{
		WeylLine[{1,p1_,m1_}, {1,p2_,m2_},WMSS[a___]]:>
			WeylLine[{"x",p1,m1},WM[a],{"xd",p2,m2}],
		WeylLine[{1,p1_,m1_}, {1,p2_,m2_}, WMSSB[a___]]:>
			WeylLine[{"x",p1,m1},WM[a],{"x",p2,m2}],
		WeylLine[{1,p1_,m1_}, {1,p2_,m2_}, WMSBS[a___]]:>
			WeylLine[{"xd",p1,m1},WM[a],{"xd",p2,m2}],
		WeylLine[{1,p1_,m1_}, {1,p2_,m2_}, WMSBSB[a___]]:>
			WeylLine[{"xd",p1,m1},WM[a],{"x",p2,m2}],

		WeylLine[{1,p1_,m1_}, {-1,p2_,m2_}, WMSS[a___]]:>
			WeylLine[{"x",p1,m1},WM[a],{"yd",p2,m2}],
		WeylLine[{1,p1_,m1_}, {-1,p2_,m2_}, WMSSB[a___]]:>
			WeylLine[{"x",p1,m1},WM[a],{"y",p2,m2}],
		WeylLine[{1,p1_,m1_}, {-1,p2_,m2_}, WMSBS[a___]]:>
			WeylLine[{"xd",p1,m1},WM[a],{"yd",p2,m2}],
		WeylLine[{1,p1_,m1_}, {-1,p2_,m2_}, WMSBSB[a___]]:>
			WeylLine[{"xd",p1,m1},WM[a],{"y",p2,m2}],

		WeylLine[{-1,p1_,m1_}, {1,p2_,m2_}, WMSS[a___]]:>
			WeylLine[{"y",p1,m1},WM[a],{"xd",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {1,p2_,m2_}, WMSSB[a___]]:>
			WeylLine[{"y",p1,m1},WM[a],{"x",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {1,p2_,m2_}, WMSBS[a___]]:>
			WeylLine[{"yd",p1,m1},WM[a],{"xd",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {1,p2_,m2_}, WMSBSB[a___]]:>
			WeylLine[{"yd",p1,m1},WM[a],{"x",p2,m2}],

		WeylLine[{-1,p1_,m1_}, {-1,p2_,m2_}, WMSS[a___]]:>
			WeylLine[{"y",p1,m1},WM[a],{"yd",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {-1,p2_,m2_}, WMSSB[a___]]:>
			WeylLine[{"y",p1,m1},WM[a],{"y",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {-1,p2_,m2_}, WMSBS[a___]]:>
			WeylLine[{"yd",p1,m1},WM[a],{"yd",p2,m2}],
		WeylLine[{-1,p1_,m1_}, {-1,p2_,m2_}, WMSBSB[a___]]:>
			WeylLine[{"yd",p1,m1},WM[a],{"y",p2,m2}]
	};

	wl=wl/.{LTS[mu__]:>LTensor[WeylS,mu]};
	wl=wl/.{LTSB[mu__]:>LTensor[WeylSBar,mu]};
	wl=wl/.{WM[a___]:>WeylMatrix[a]};
	wl
];


InternalFermionSpinSum[WeylLine[l1___],WeylLine[l2___]]:=Module[{wl1,wl2,prod,WFP,LD},

	LD[a_,b_]:=LDot[a,b];
	wl1=ConvertWeylLine[WeylLine[l1]];
	wl2=ConvertWeylLine[WeylLine[l2]];

	prod={wl1,wl2}/.{
		{WeylLine[{s1_,p1_,m1_},WeylMatrix[a___],{s2_,p2_,m2_}],
		 WeylLine[{s3_,p2_,m2_},WeylMatrix[b___],{s4_,p1_,m1_}]} :>
		WeylMatrix[a,WFP[s2<>s3,p2,m2],b,WFP[s4<>s1,p1,m1]]
	};

	prod=prod//.{
		WeylMatrix[a___,WFP["xxd",p_,m_],b___]:>WeylMatrix[a,LD[WeylS,p],b],
		WeylMatrix[a___,WFP["xdx",p_,m_],b___]:>WeylMatrix[a,LD[WeylSBar,p],b],
		WeylMatrix[a___,WFP["yyd",p_,m_],b___]:>WeylMatrix[a,LD[WeylS,p],b],
		WeylMatrix[a___,WFP["ydy",p_,m_],b___]:>WeylMatrix[a,LD[WeylSBar,p],b],
		WeylMatrix[a___,WFP["xy",p_,m_],b___]:>WeylMatrix[a,m*Weyl1,b],
		WeylMatrix[a___,WFP["xdyd",p_,m_],b___]:>WeylMatrix[a,-m*Weyl1,b],
		WeylMatrix[a___,WFP["yx",p_,m_],b___]:>WeylMatrix[a,-m*Weyl1,b],
		WeylMatrix[a___,WFP["ydxd",p_,m_],b___]:>WeylMatrix[a,m*Weyl1,b]
	};

	WeylTrace[prod]
];


FermionSpinSum[exp_]:=Module[{iexp,lines,nlines,coeffs},

	ExpandAll[exp]/.{
		WeylLine[
			wf1_List /; SameQ[Length[wf1], 3],
			wf2_List /; SameQ[Length[wf2], 3], WeylMatrix[mtx1___]
		] *
		WeylLine[
			wf3_List /; SameQ[Length[wf3], 3],
			wf4_List /; SameQ[Length[wf4], 3], WeylMatrix[mtx2___]
		] :> InternalFermionSpinSum[WeylLine[wf1,wf2,WeylMatrix[mtx1]],
							        WeylLine[wf3,wf4,WeylMatrix[mtx2]]]
	}
]

End[]
