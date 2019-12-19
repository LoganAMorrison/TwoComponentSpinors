(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)


WeylTrace::usage =
 "WeylTrace[WeylMatrix[\!\(\*SubscriptBox[\(mtx\), \(1\)]\),...,\!\(\
\*SubscriptBox[\(mtx\), \(n\)]\)]] computes the trace of the product \
of sigma matrices \!\(\*SubscriptBox[\(mxt\), \(1\)]\)\[CenterDot]\
\[CenterDot]\[CenterDot]\!\(\*SubscriptBox[\(mxt\), \(n\)]\).";

WeylTrace::InvalidOddNumWeylArgs="Odd number of WeylS and WeylSBar
encountered in WeylTrace[...].";


(* Repeated WeylS or WeylSBar *)
WeylTrace::InvalidSpinorStructureSS="Invalid spinor structure of \ WelyMatrix[...,WeylS,WeylS,...] encountered in WeylTrace[...].";

WeylTrace::InvalidSpinorStructureWelySBSB="Invalid spinor structure \ of WelyMatrix[...,WeylSBar, WeylSBar,...] encountered in WeylTrace[...].";


(* WeylSBar * Weyl1Bar * ... or WeylS * Weyl1 * ... *)
WeylTrace::InvalidSpinorStructureSB1B="Invalid spinor structure \ of WelyMatrix[...,WeylSBar,Weyl1Bar,...] encountered in WeylTrace[...].";

WeylTrace::InvalidSpinorStructureS1="Invalid spinor structure \ of WelyMatrix[...,WeylS,Weyl1,...] encountered in WeylTrace[...].";


(* repeated Weyl1 or Weyl1Bar*)
WeylTrace::InvalidSpinorStructure11B="Invalid spinor structure \ of WelyMatrix[...,Weyl1,Weyl1Bar,...] encountered in WeylTrace[...].";

WeylTrace::InvalidSpinorStructure1B1="Invalid spinor structure \ of WelyMatrix[...,Weyl1Bar,Weyl1,...] encountered in WeylTrace[...].";

Begin["Private`"]


LT[a__]:=LTensor[a];
WM[a___]:=WeylMatrix[a];
MT=MetricG;

(*
    ValidateWeylMatrix[WeylMatrix[mtx1,...,mtxn]]

Make sure the arguments of 'WeylMatrix' are valid. Valid arguments must be
alternating in WeylS and WeylSBar. Combinations of
    WeylMatrix[..., WeylSBar, Weyl1, WeylS,...]
    WeylMatrix[..., WeylS, Weyl1Bar, WeylSBar,...]
are acceptable. Combinations of 'WeylSBar, Weyl1Bar, WeylS' or
'WeylS, Weyl1, WeylSBar' are NOT acceptable.
*)
ValidateWeylMatrix[WeylMatrix[args___]]:=Module[{expr,iargs},

    expr=WeylMatrix[args];

    (* check for WeylMatrix[...,WeylS, WeylS,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,iLTS[_],iLTS[_],right___]],
        Message[WeylTrace::InvalidSpinorStructureSS];
        Print[WeylMatrix@@args];Abort[];
    ];
    (* check for WeylMatrix[...,WeylSBar, WeylSBar,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,iLTSB[_],iLTSB[_],right___]],
        Message[WeylTrace::InvalidSpinorStructureWelySBSB];
        Print[WeylMatrix[args]];Abort[];
    ];

    (* check for WeylMatrix[...,WeylSBar, Weyl1Bar,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,iLTSB[_],Weyl1Bar,right___]],
        Message[WeylTrace::InvalidSpinorStructureSB1B];
        Print[WeylMatrix[args]];Abort[];
    ];

    (* check for WeylMatrix[...,WeylS, Weyl1, WeylSBar,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,iLTS[_],Weyl1,right___]],
        Message[WeylTrace::InvalidSpinorStructureS1];
        Print[WeylMatrix[args]];Abort[];
    ];

    (* check for WeylMatrix[...,Weyl1, Weyl1Bar,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,Weyl1,Weyl1Bar,right___]],
        Message[WeylTrace::InvalidSpinorStructure11B];
        Print[WeylMatrix[args]];Abort[];
    ];
    (* check for WeylMatrix[...,Weyl1Bar, Weyl1,...]*)
    If[MatchQ[expr,
        WeylMatrix[left___,Weyl1Bar,Weyl1,right___]],
        Message[WeylTrace::InvalidSpinorStructure1B1];
        Print[WeylMatrix[args]];Abort[];
    ];

    iargs=DeleteCases[List[args], Weyl1|Weyl1Bar];
    (* throw if there are an odd number of pauli matrices *)
    If[OddQ[Length[iargs]],
        Message[WeylTrace::InvalidOddNumWeylArgs];
        Print[WeylMatrix[args]];Abort[];
    ];

    Return["Valid"];
];

WeylTrace[WeylMatrix[args___]]:=Module[{exp,clist,wlist,i},

	exp=ExpandAll[WM[args]];
	(* uncontract sigma matrices *)
	exp=X`Utilities`Uncontract[#, WeylS|WeylSBar]&/@exp;

	(* relabel LTensor[WeylS, mu] and LTensor[WeylSBar, mu] *)
	exp=exp/.{LTensor[WeylS, mu_]:>iLTS[mu],LTensor[WeylSBar, mu_]:>iLTSB[mu]};
    (* split internal sums: Weyl[a+b,...] -> Weyl[a,...] + Weyl[b,...] *)
    exp=Distribute[exp, Plus, WeylMatrix];
    (* pull out constants *)
    exp=exp//.{
        WM[left___, a__ * iLTS[mu_], right___] :> a * WM[left, iLTS[mu], right],
        WM[left___, a__ * iLTSB[mu_], right___] :> a * WM[left, iLTSB[mu], right],
        WM[left___, a__ * Weyl1, right___] :> a * WM[left, Weyl1, right],
        WM[left___, a__ * Weyl1Bar, right___] :> a * WM[left, Weyl1Bar, right]
    };

	(* put
		a*WM1[mtx11,...,mtx1n] + b*WM2[mtx21,...,mtx2n] + ...
	into
		clist = {a,b,...}
	and
		wlist = {{mtx11,...,mtx1n},{mtx21,...,mtx2n},...}
	*)
	exp=If[SameQ[Head[exp],Plus], List@@exp, {exp}];
	clist=Table[elm/.{c_. * WM[b___]:> c}, {elm, exp}];
	wlist=Table[elm/.{c_. * WM[b___]:> WM[b]}, {elm, exp}];

    (* go through wlist and make sure all is valid. Will abort if invalid *)
    Map[ValidateWeylMatrix,wlist];

    (* Remove Weyl1's and Weyl1Bar where not needed *)
    wlist=wlist//.{
        (* expressions like WM[L,a*1] or WM[L,a*1,R] *)
        WM[left__, a_. * Weyl1, right___] :> a * WM[left, right],
        (* expressions like WM[a*1,R] or WM[L,a*1,R] *)
        WM[left___, a_. * Weyl1, right__] :> a * WM[left, right],
        (* expressions like WM[L,a*1b] or WM[L,a*1b,R] *)
        WM[left__, a_. * Weyl1Bar, right___] :> a * WM[left, right],
        (* expressions like WM[a*1b] or WM[R,a*1b,L] *)
        WM[left___, a_. * Weyl1Bar, right__] :> a * WM[left, right]
    };

	(* Evaluate traces by recursion *)
    wlist = ReplaceRepeated[wlist,{
        (* base case of Tr[1] *)
        WM[] :> 2,
        WM[ Weyl1 ] :> 2,
        WM[ Weyl1Bar ] :> 2,
		(* base case of Tr[s, sb] *)
		WM[ iLTS[mu1_], iLTSB[mu2_] ] :> 2 * LT[MT, mu1, mu2],
		(* base case of Tr[sb, b] *)
		WM[ iLTSB[mu1_], iLTS[mu2_] ] :> 2 * LT[MT, mu1, mu2],
		(* base case of Tr[sb, s, sb, s] *)
		WM[ iLTS[mu1_], iLTSB[mu2_], iLTS[mu3_], iLTSB[mu4_] ] :> 2 * (
			LT[MT, mu1, mu2]*LT[MT, mu3, mu4] -
			LT[MT, mu1, mu3]*LT[MT, mu2, mu4] +
			LT[MT, mu1, mu4]*LT[MT, mu2, mu3] +
			I * LT[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* base case of Tr[s, sb, s, sb] *)
		WM[ iLTSB[mu1_], iLTS[mu2_], iLTSB[mu3_], iLTS[mu4_] ] :> 2 * (
			LT[MT, mu1, mu2]*LT[MT, mu3, mu4] -
			LT[MT, mu1, mu3]*LT[MT, mu2, mu4] +
			LT[MT, mu1, mu4]*LT[MT, mu2, mu3] -
			I * LT[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* reduce Tr[s, sb, s, ....] to sum of Tr[s, ....] *)
		WM[ iLTS[mu1_], iLTSB[mu2_], iLTS[mu3_], a__ ] :> (
			LT[MT, mu1, mu2] * WM[ iLTS[mu3], a ] -
			LT[MT, mu1, mu3] * WM[ iLTS[mu2], a ] +
			LT[MT, mu2, mu3] * WM[ iLTS[mu1], a ] +
			I * LT[X`LeviCivitaE, mu1, mu2, mu3, mu] * WM[ iLTS[mu], a ]
		)/.{mu->Unique["$"]},
		(* reduce Tr[sb, s, sb, ....] to sum of Tr[sb, ....] *)
		WM[ iLTSB[mu1_], iLTS[mu2_], iLTSB[mu3_], a__ ] :> (
			LT[MT, mu1, mu2] * WM[ iLTSB[mu3], a ] -
			LT[MT, mu1, mu3] * WM[ iLTSB[mu2], a ] +
			LT[MT, mu2, mu3] * WM[ iLTSB[mu1], a ] -
			I * LT[X`LeviCivitaE, mu1, mu2, mu3, mu] * WM[ iLTSB[mu], a ]
		)/.{mu->Unique["$"]}
	}];

	(* put coefficients back in *)
	Dot[clist, wlist]
];


WeylTrace[a___] := WeylTrace[WeylMatrix[a]];


End[]
