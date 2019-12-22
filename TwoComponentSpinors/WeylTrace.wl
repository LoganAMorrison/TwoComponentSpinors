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

WeylTraceInternal[(WM:iWMLL|iWMRR)[args___]]:=Module[{mu},
    mu=Unique["mu$"];
    ReplaceAll[WM[args],{
        (* base case of Tr[1] *)
        WM[] :> 2,
		(* base case of Tr[s, sb] or Tr[sb, s] *)
		WM[ iLTS[mu1_], iLTS[mu2_] ] :> 2 * LTensor[MetricG, mu1, mu2],
		(* base case of Tr[sb, s, sb, s] or Tr[s, sb, s, sb] *)
		WM[ iLTS[mu1_], iLTS[mu2_], iLTS[mu3_], iLTS[mu4_] ] :>
        2 * (
			LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu4] -
			LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu4] +
			LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu3] +
            Switch[WM,iWMLL,1,iWMRR,-1] *
			I * LTensor[LeviCivitaE, mu1, mu2, mu3, mu4]
		),
		(* reduce Tr[X, s, sb, s] or Tr[X, sb, s, sb] to sum of Tr[X, s(b)] *)
		WM[left__, iLTS[mu1_], iLTS[mu2_], iLTS[mu3_]] :>
        (
			LTensor[MetricG, mu1, mu2] * WM[left, iLTS[mu3]] -
			LTensor[MetricG, mu1, mu3] * WM[left, iLTS[mu2]] +
			LTensor[MetricG, mu2, mu3] * WM[left, iLTS[mu1]] +
            Switch[WM,iWMLL,-1,iWMRR,1] *
			I * LTensor[LeviCivitaE, mu1, mu2, mu3, mu] * WM[left, iLTS[mu]]
		)
    }]
];

WeylTrace[(WM:WeylMatrixL|WeylMatrixR)[args___]]:=Module[{iexpr},
    iexpr=ConvertToInternal[WeylMatrixExpand[WM[args]]];
    (* make sure there are an even number of sigma matrices, in which case
       the head will be iWMRR or iWMLL. *)
    If[Or[SameQ[Head[iexpr], iWMRL],SameQ[Head[iexpr], iWMLR]],
        Message[WeylTrace::InvalidOddNumWeylArgs];Abort[];
    ];
    (* apply WeylTraceInternal to all iWMLL and iWMRR *)
    iexpr=iexpr//.{
        mtx_iWMLL :> WeylTraceInternal[mtx],
        mtx_iWMRR :> WeylTraceInternal[mtx]
    };
    (* put coefficients back in *)
	ConvertToExternal[iexpr]
];

(* CODE THAT SHOULD BE TAKEN CARE OF BY WeylMatrixExpand

uncontract sigma matrices
iexpr=ConvertToInternal[X`Utilities`Uncontract[#, WeylS]&/@iexpr];

split internal sums: Weyl[a+b,...] -> Weyl[a,...] + Weyl[b,...]
iexpr=Distribute[iexpr, Plus, iWMLL];
iexpr=Distribute[iexpr, Plus, iWMRR];

pull out constants
iexpr=iexpr//.{
    (pat:iWMLL|iWMRR)[left___, a__ * iLTS[mu_], right___] :>
        a * pat[left, iLTS[mu], right]
};
*)

End[]
