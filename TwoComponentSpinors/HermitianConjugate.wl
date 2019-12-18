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

Options[HermitianConjugate] = {
	"ComplexSymbols" -> {}
}

HermitianConjugate[exp_, OptionsPattern[]]:= Module[{iexp,cvars},

	cvars = OptionValue["ComplexSymbols"];
	(* Validate that "ComplexSymbols" is a list *)
	If[Not[ListQ[cvars]], $Failed];
	(* make replacement list for complex valiables *)
	cvars = Table[cvar -> Conjugate[cvar], {cvar, cvars}];

	iexp=exp/.{I->-I};
	iexp=iexp/.{
		PolarizationVector[p_,m_]:>PolarizationVectorDag[p,m],
		PolarizationVectorDag[p_,m_]:>PolarizationVector[p,m]
	};
	iexp=iexp/.{WeylMatrix[a___]:>WeylMatrix@@Reverse[List[a]]};
	iexp=iexp/.{WeylLine[a_List,b_List,c_WeylMatrix]:>WeylLine[b,a,c]};
	iexp/.cvars
]

End[]
