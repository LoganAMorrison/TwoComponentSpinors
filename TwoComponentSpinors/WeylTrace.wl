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

WeylTrace::InvalidOddNumWeylArgs="Odd number of sigma matrices encountered in\
WeylTrace[...] or perhaps a missing Weyl1.";

WeylTrace::InvalidSpinorStructure="Invalid spinor structure encountered in \
WeylTrace[...]. Either two adjacent WeylS's, two adjacent WeylSBar's or \
non-wely matrix encountered.";

Begin["Private`"]


WeylTrace[WeylMatrix[args___]]:=Module[{exp,WM,LTS,LTSB,LT,MT,clist,wlist,i,iargs},

	LT[a__]:=LTensor[a];
    WM[a___]:=WeylMatrix[a];
	MT=MetricG;

	exp=ExpandAll[WM[args]];
	(* uncontract sigma matrices *)
	exp=X`Utilities`Uncontract[exp, WeylS];
	exp=X`Utilities`Uncontract[exp, WeylSBar];

	(* relabel LTensor[WeylS, mu] and LTensor[WeylSBar, mu]*)
	exp=exp/.{LTensor[WeylS, mu_]:>LTS[mu]};
	exp=exp/.{LTensor[WeylSBar, mu_]:>LTSB[mu]};

	(* replace inner WeylMatrix[args] with args *)
	exp=exp//.{WM[a___,WM[b___],c___]:>WM[a,b,c]};

	(* split sums *)
    exp=exp//.{
        WM[mtx1___, a_ + b_, mtx2___] :> WM[mtx1, a, mtx2] + WM[mtx1, b, mtx2]
    };
    (* pull out constants *)
    exp=exp//.{
        WM[mtx1___, a__ * LTS[mu1_], mtx2___] :> a * WM[mtx1, LTS[mu1], mtx2],
        WM[mtx1___, a__ * LTSB[mu1_], mtx2___] :> a * WM[mtx1, LTSB[mu1], mtx2]
    };
    (* Remove Weyl1's *)
    exp=exp//.{
        WM[mtx1___, a_. * Weyl1, mtx2___]:>a*WM[mtx1,mtx2]
    };

	exp=ExpandAll[exp];

	(** VALIDATION **)

	(* put
		a*WM1[mtx11,...,mtx1n] + b*WM2[mtx21,...,mtx2n] + ...
	into
		clist = {a,b,...}
	and
		wlist = {{mtx11,...,mtx1n},{mtx21,...,mtx2n},...}
	*)
	exp=If[SameQ[Head[exp],Plus], List@@exp, {exp}];
	clist=Table[elm/.{c_. * WM[b___]:> c}, {elm, exp}];
	wlist=Table[elm/.{c_. * WM[b___]:> List[b]}, {elm, exp}];\

	(* check that arguments are alternating *)
	For[i=1, i<=Length[wlist],i++,
		iargs=DeleteCases[DeleteCases[wlist[[i]],Weyl1],Weyl1Bar];
		(* throw if there are an odd number of pauli matrices *)
		If[OddQ[Length[iargs]],
			Message[WeylTrace::InvalidOddNumWeylArgs];Return[$Failed];
		];

		evenargs=DeleteDuplicates[
			Table[Head[iargs[[2*i]]],{i,Floor[Length[iargs]/2]}]
		];
		oddargs=DeleteDuplicates[
			Table[Head[iargs[[2*i-1]]],{i,Ceiling[Length[iargs]/2]}]
		];

        Print[evenargs];
        Print[oddargs];

		(* lengths of evenargs and oddargs MUST be 1 and evenargs!=oddargs *)
		If[Or[
			And[UnsameQ[Length[evenargs],1], UnsameQ[Length[evenargs],0]],
			And[UnsameQ[Length[evenargs],1], UnsameQ[Length[evenargs],0]],
			SameQ[evenargs[[1]],oddargs[[1]]]
		],
			Message[WeylTrace::InvalidSpinorStructure];Return[$Failed];
		];

	];

	(* put back into WM *)
	wlist=Table[WM@@iargs, {iargs, wlist}];

    wlist/.{WM[a___, Weyl1, b___]:> WM[a,b]};
    wlist/.{WM[a___, Weyl1Bar, b___]:> WM[a,b]};

	(* otherwise, recurrsively evaluate *)
	wlist = ReplaceRepeated[wlist,{
        (* base case of Tr[1] *)
        WM[ ] :> 2,
        WM[ Weyl1 ] :> 2,
        WM[ Weyl1Bar ] :> 2,
		(* base case of Tr[s, sb] *)
		WM[ LTS[mu1_], LTSB[mu2_] ] :> 2 * LT[MT, mu1, mu2],
		(* base case of Tr[sb, b] *)
		WM[ LTSB[mu1_], LTS[mu2_] ] :> 2 * LT[MT, mu1, mu2],
		(* base case of Tr[sb, s, sb, s] *)
		WM[ LTS[mu1_], LTSB[mu2_], LTS[mu3_], LTSB[mu4_] ] :> 2 * (
			LT[MT, mu1, mu2]*LT[MT, mu3, mu4] -
			LT[MT, mu1, mu3]*LT[MT, mu2, mu4] +
			LT[MT, mu1, mu4]*LT[MT, mu2, mu3] +
			I * LT[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* base case of Tr[s, sb, s, sb] *)
		WM[ LTSB[mu1_], LTS[mu2_], LTSB[mu3_], LTS[mu4_] ] :> 2 * (
			LT[MT, mu1, mu2]*LT[MT, mu3, mu4] -
			LT[MT, mu1, mu3]*LT[MT, mu2, mu4] +
			LT[MT, mu1, mu4]*LT[MT, mu2, mu3] -
			I * LT[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* reduce Tr[s, sb, s, ....] to sum of Tr[s, ....] *)
		WM[ LTS[mu1_], LTSB[mu2_], LTS[mu3_], a__ ] :> (
			LT[MT, mu1, mu2] * WM[ LTS[mu3], a ] -
			LT[MT, mu1, mu3] * WM[ LTS[mu2], a ] +
			LT[MT, mu2, mu3] * WM[ LTS[mu1], a ] +
			I * LT[X`LeviCivitaE, mu1, mu2, mu3, mu] * WM[ LTS[mu], a ]
		)/.{mu->Unique["$"]},
		(* reduce Tr[sb, s, sb, ....] to sum of Tr[sb, ....] *)
		WM[ LTSB[mu1_], LTS[mu2_], LTSB[mu3_], a__ ] :> (
			LT[MT, mu1, mu2] * WM[ LTSB[mu3], a ] -
			LT[MT, mu1, mu3] * WM[ LTSB[mu2], a ] +
			LT[MT, mu2, mu3] * WM[ LTSB[mu1], a ] -
			I * LT[X`LeviCivitaE, mu1, mu2, mu3, mu] * WM[ LTSB[mu], a ]
		)/.{mu->Unique["$"]}
	}];

	(* put coefficients back in *)
	Dot[clist, wlist]
];


WeylTrace[a___] := WeylTrace[WeylMatrix[a]];


End[]
