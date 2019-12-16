(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

HermitianConjugate::usage="HermitianConjugate[x] returns the hermitian conjugate of x.";


Begin["Private`"]

HermitianConjugate[exp_]:= ReplaceAll[Expand[exp],{
	(* polarization vectors *)
	PolarizationVector[p_,m_]:>PolarizationVectorDag[p,m],
	PolarizationVectorDag[p_,m_]:>PolarizationVector[p,m],
	(* spinor wave functions *)
	SpinorX[p_,m_]:>SpinorXDag[p,m],
	SpinorXDag[p_,m_]:>SpinorX[p,m],
	SpinorY[p_,m_]:>SpinorYDag[p,m],
	SpinorYDag[p_,m_]:>SpinorY[p,m],
	SigmaMatrix[a__]:>SigmaMatrix@@Reverse[List[a]],
	SpinorLine[a_, SigmaMatrix[b__], c_] :> SpinorLine[
		HermitianConjugate[c],
		HermitianConjugate[SigmaMatrix[b]],
		HermitianConjugate[a]
	]
}];

End[]

EndPackage[]
