(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

PolarizationSum::usage="PolarizationSum[exp] Perform the sum over \
polarization vectors.";
Bar::usage="Bar[k] If k = (E, kk), then Bar[k] = (E, -kk) where kk is the \
vector component of k.";

Begin["Private`"]

PolarizationSum[exp_]:=Module[{uncontracted},
    (* Continure to uncontract until all polarization vectors are bare *)
	uncontracted=FixedPoint[UncontractPolarizationVector,exp];

	ReplaceAll[uncontracted,{
		X`LTensor[PolarizationVector[k_,0],mu1_]*X`LTensor[PolarizationVectorDag[k_,0],mu2_]:>-X`LTensor[X`MetricG,mu1,mu2] + (X`LTensor[k,mu1]*X`LTensor[Bar[k],mu2] + X`LTensor[k,mu2]*X`LTensor[Bar[k],mu1]) / X`LDot[k,Bar[k]],
		X`LTensor[PolarizationVector[k_,m_],mu1_]*X`LTensor[PolarizationVectorDag[k_,m_],mu2_]:>-X`LTensor[X`MetricG,mu1,mu2]+X`LTensor[k,mu1]X`LTensor[k,mu2]/m^2
	}]
]

End[]

EndPackage[]
