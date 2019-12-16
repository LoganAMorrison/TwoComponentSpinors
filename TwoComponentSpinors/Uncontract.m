(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

Begin["Private`"]

UncontractSigmaMatrix[exp_]:=Module[{mu},
	exp/.{
	SigmaMatrix[a___,X`LDot[x_,PauliS],b___]:>X`LTensor[x,mu]*SigmaMatrix[a,X`LTensor[PauliS,mu],b],
	SigmaMatrix[a___,X`LDot[PauliS, x_],b___]:>X`LTensor[x,mu]*SigmaMatrix[a,X`LTensor[PauliS,mu],b],
	SigmaMatrix[a___,X`LDot[x_,PauliSBar],b___]:>X`LTensor[x,mu]*SigmaMatrix[a,X`LTensor[PauliSBar,mu],b],
	SigmaMatrix[a___,X`LDot[PauliSBar,x_],b___]:>X`LTensor[x,mu]*SigmaMatrix[a,X`LTensor[PauliSBar,mu],b]
	}
];

UncontractPolarizationVector[exp_]:=Module[{mu},
	exp/.{
		X`LDot[PolarizationVector[k_,m_],x_]:>X`LTensor[x,mu]*X`LTensor[PolarizationVector[k,m],mu],
		X`LDot[x_,PolarizationVector[k_,m_]]:>X`LTensor[x,mu]*X`LTensor[PolarizationVector[k,m],mu],
		X`LDot[PolarizationVectorDag[k_,m_],x_]:>X`LTensor[x,mu]*X`LTensor[PolarizationVectorDag[k,m],mu],
		X`LDot[x_,PolarizationVectorDag[k_,m_]]:>X`LTensor[x,mu]*X`LTensor[PolarizationVectorDag[k,m],mu]
	}
];

End[]

EndPackage[]
