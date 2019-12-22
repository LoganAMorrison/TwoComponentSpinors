(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

HermitianConjugate::usage="HermitianConjugate[x] returns the hermitian conjugate of x.";

HermitianConjugate::InvalidComplexSymbols="ComplexSymbols must be a list of \
Symbol's.";


Begin["Private`"]

HermitianConjugateWeylMatrix[(WM:WeylMatrixR|WeylMatrixL)[args___]]:=Module[{iexpr},
	iexpr=ConvertToInternal[WeylMatrixExpand[WM[args]]];
	(* check if expression is RR *)
	If[Not[FreeQ[iexpr,iWMRR]],
		Return[ConvertToExternal[iWMLL@@Reverse[List[args]]]]
	];
	(* check if expression is LL *)
	If[Not[FreeQ[iexpr,iWMLL]],
		Return[ConvertToExternal[iWMRR@@Reverse[List[args]]]]
	];
	(* check if expression is RL *)
	If[Not[FreeQ[iexpr,iWMRL]],
		Return[ConvertToExternal[iWMRL@@Reverse[List[args]]]]
	];
	(* check if expression is LR *)
	If[Not[FreeQ[iexpr,iWMLR]],
		Return[ConvertToExternal[iWMLR@@Reverse[List[args]]]]
	];
]

Options[HermitianConjugate] = {
	"ComplexSymbols" -> {}
}

HermitianConjugate[expr_, OptionsPattern[]]:= Module[{iexpr,cvars},

	cvars = OptionValue["ComplexSymbols"];
	(* Validate that "ComplexSymbols" is a list *)
	If[Not[ListQ[cvars]], $Failed];
	(* make replacement list for complex valiables *)
	cvars=Table[cvar -> Conjugate[cvar], {cvar, cvars}];

	(* Change a+ib -> a-ib*)
	iexpr=expr/.{Complex[a_,b_]:>Complex[a,-b]};

	iexpr=iexpr/.{
		PolarizationVector[p_,m_]:>Conjugate[PolarizationVector[p,m]]
	};
	iexpr=iexpr/.{
		(pat:WeylMatrixL|WeylMatrixR)[a___] :>
			HermitianConjugateWeylMatrix[pat[a]]
	};
	iexpr=iexpr/.{
		WeylLine[a_List,b_List,(pat:WeylMatrixL|WeylMatrixR)[mtx___]] :>
			WeylLine[b,a,pat[mtx]]
	};
	iexpr/.cvars
]

End[]
