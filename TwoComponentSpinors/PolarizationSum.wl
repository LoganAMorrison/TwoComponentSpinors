(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

PolarizationSum::usage =
  "PolarizationSum[exp, k] or PolarizationSum[exp, \
{\!\(\*SubscriptBox[\(k\), \(1\)]\),...,\!\(\*SubscriptBox[\(k\), \(n\
\)]\)}] performs the sum over polarizations for polarization vectors \
k or \!\(\*SubscriptBox[\(k\), \(1\)]\),...,\!\(\*SubscriptBox[\(k\), \
\(n\)]\).";

(*
Bar::usage="Bar[k] If k = (E, kk), then Bar[k] = (E, -kk) where kk is the \
vector component of k.";
*)

Begin["`Private`"]

PolarizationSum[expr_, k_]:=Module[{iexpr},

	iexpr=X`Utilities`Uncontract[expr, PolarizationVector[k,_]];

	ReplaceAll[iexpr,{
		(* massless polarization spin sum *)
		LTensor[PolarizationVector[k,0],mu_] *
		LTensor[Conjugate[PolarizationVector[k,0]],nu_] :>
			-LTensor[MetricG,mu,nu],

		(* massive polarization spin sum *)
		LTensor[PolarizationVector[k,m_],mu_] *
		LTensor[Conjugate[PolarizationVector[k,m_]],nu_]:>
			-LTensor[MetricG,mu,nu] + LTensor[k,mu]*LTensor[k,nu] / m^2
	}]
]

PolarizationSum[expr_,ks_List]:=Module[{iexpr,i},
	iexpr=Expand[expr];
	(* uncontract all momenta in ks*)
	For[i=1,i<=Length[ks],i++,iexpr=X`Utilities`Uncontract[iexpr, PolarizationVector[ks[[i]]]]];

	iexpr//.{
		(* massless polarization spin sum *)
		LTensor[PolarizationVector[k_,0],mu_] *
		LTensor[Conjugate[PolarizationVector[k_,0]],nu_] :>
			-LTensor[MetricG,mu,nu]/;MemberQ[ks,k],

		(* massive polarization spin sum *)
		LTensor[PolarizationVector[k_,m_],mu_] *
		LTensor[Conjugate[PolarizationVector[k_,m_]],nu_]:>
			-LTensor[MetricG,mu,nu] +
			LTensor[k,mu]*LTensor[k,nu] / m^2/;MemberQ[ks,k]
	}
]

End[]
