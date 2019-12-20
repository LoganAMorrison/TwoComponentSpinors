(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)


WeylTrace::usage="WeylTrace[WeylMatrixL[\!\(\*SubscriptBox[\(mtx\), \
\(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)]] or \
WeylTrace[WeylMatrixR[\!\(\*SubscriptBox[\(mtx\), \
\(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)]] computes the \
trace of the product of sigma matrices (\!\(\*SubscriptBox[\(mtx\), \
\(1\)]\)\[Times]\[CenterDot]\[CenterDot]\[CenterDot]\[Times]\!\(\*\
SubscriptBox[\(mtx\), \(n\)]\)).";

WeylTrace::InvalidOddNumWeylArgs="Odd number of Sigma matrices \
encountered in WeylTrace[...].";

Begin["Private`"]


LT[a__]:=LTensor[a];
MT=MetricG;

(*
    ValidateWeylMatrix[WeylMatrix[mtx1,...,mtxn]]

Make sure the arguments of 'WeylMatrixL(R)' are valid. There must be an even
number of WeylS's.
*)
ValidateWeylMatrix[(pat:WeylMatrixL|WeylMatrixR)[args___]]:=Module[{expr,iargs},

    expr=pat[args];

    iargs=DeleteCases[List[args], Weyl1];
    (* throw if there are an odd number of pauli matrices *)
    If[OddQ[Length[iargs]],
        Message[WeylTrace::InvalidOddNumWeylArgs];
        Print[pat[args]];Abort[];
    ];

    Return["Valid"];
];

WeylTrace[(WM:WeylMatrixL|WeylMatrixR)[args___]]:=Module[{exp,clist,wlist,i},

	exp=ExpandAll[ConvertToInternal[WM[args]]];
	(* uncontract sigma matrices *)
	exp=X`Utilities`Uncontract[#, WeylS]&/@exp;

	(* relabel LTensor[WeylS, mu] *)
	exp=exp/.{LTensor[WeylS, mu_]:>iLTS[mu]};
    (* split internal sums: Weyl[a+b,...] -> Weyl[a,...] + Weyl[b,...] *)
    exp=Distribute[exp, Plus, iWML];
    exp=Distribute[exp, Plus, iWMR];
    (* pull out constants *)
    exp=exp//.{
        (pat:iWML|iWMR)[left___, a__ * iLTS[mu_], right___] :>
            a * pat[left, iLTS[mu], right],
        (pat:iWML|iWMR)[left___, a__ * Weyl1, right___] :>
            a * pat[left, Weyl1, right]
    };

	(* put
		a*WM1[mtx11,...,mtx1n] + b*WM2[mtx21,...,mtx2n] + ...
	into
		clist = {a,b,...}
	and
		wlist = {{mtx11,...,mtx1n},{mtx21,...,mtx2n},...}
	*)
	exp=If[SameQ[Head[exp],Plus], List@@exp, {exp}];
	clist=Table[elm/.{c_. * (pat:iWML|iWMR)[mtx___]:> c}, {elm, exp}];
	wlist=Table[elm/.{c_. * (pat:iWML|iWMR)[mtx___]:> pat[mtx]}, {elm, exp}];

    (* go through wlist and make sure all is valid. Will abort if invalid *)
    Map[ValidateWeylMatrix,wlist];

    (* Remove Weyl1's and Weyl1Bar where not needed *)
    wlist=wlist//.{
        (* expressions like WM[L,a*1] or WM[L,a*1,R] *)
        (pat:iWML|iWMR)[left__, a_. * Weyl1, right___] :> a * pat[left, right],
        (* expressions like WM[a*1,R] or WM[L,a*1,R] *)
        (pat:iWML|iWMR)[left___, a_. * Weyl1, right__] :> a * pat[left, right]
    };

	(* Evaluate traces by recursion *)
    wlist = ReleaseHold[ReplaceRepeated[wlist,{
        (* base case of Tr[1] *)
        (pat:iWML|iWMR)[] :> 2,
        (pat:iWML|iWMR)[ Weyl1 ] :> 2,
		(* base case of Tr[s, sb] or Tr[sb, s] *)
		(pat:iWML|iWMR)[ iLTS[mu1_], iLTS[mu2_] ] :> 2 * LT[MT, mu1, mu2],
		(* base case of Tr[sb, s, sb, s] or Tr[s, sb, s, sb] *)
		(pat:iWML|iWMR)[ iLTS[mu1_], iLTS[mu2_], iLTS[mu3_], iLTS[mu4_] ] :>
        2 * (
			LT[MT, mu1, mu2]*LT[MT, mu3, mu4] -
			LT[MT, mu1, mu3]*LT[MT, mu2, mu4] +
			LT[MT, mu1, mu4]*LT[MT, mu2, mu3] +
            Hold[Switch][pat,iWML,1,iWMR,-1] *
			I * LT[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* reduce Tr[X, s, sb, s] or Tr[X, sb, s, sb] to sum of Tr[X, s(b)] *)
		(pat:iWML|iWMR)[left__, iLTS[mu1_], iLTS[mu2_], iLTS[mu3_]] :>
        (
			LT[MT, mu1, mu2] * pat[left, iLTS[mu3]] -
			LT[MT, mu1, mu3] * pat[left, iLTS[mu2]] +
			LT[MT, mu2, mu3] * pat[left, iLTS[mu1]] +
            Hold[Switch][pat,iWML,-1,iWMR,1] *
			I * LT[X`LeviCivitaE, mu1, mu2, mu3, mu] * pat[left, iLTS[mu]]
		)/.{mu->Unique["$"]}
	}]];

	(* put coefficients back in *)
	ConvertToExternal[Dot[clist, wlist]]
];

End[]
