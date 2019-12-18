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

HermitianConjugate[expr_, OptionsPattern[]]:= Module[{iexpr,cvars},

	cvars = OptionValue["ComplexSymbols"];
	(* Validate that "ComplexSymbols" is a list *)
	If[Not[ListQ[cvars]], $Failed];
	(* make replacement list for complex valiables *)
	cvars=Join[{I->-I,-I->I},Table[cvar -> Conjugate[cvar], {cvar, cvars}]];

	iexpr=expr/.{Weyl1->Weyl1Bar, Weyl1Bar->Weyl1};
	iexpr=iexpr/.{
		PolarizationVector[p_,m_]:>PolarizationVectorDag[p,m],
		PolarizationVectorDag[p_,m_]:>PolarizationVector[p,m]
	};
	iexpr=iexpr/.{WeylMatrix[a___]:>WeylMatrix@@Reverse[List[a]]};
	iexpr=iexpr/.{WeylLine[a_List,b_List,c_WeylMatrix]:>WeylLine[b,a,c]};
	iexpr/.cvars
]

End[]
