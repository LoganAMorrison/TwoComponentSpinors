(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

PolarizationSum::usage="PolarizationSum[exp, k] Perform the sum over \
polarization vectors corresponding to momentum k.";

Bar::usage="Bar[k] If k = (E, kk), then Bar[k] = (E, -kk) where kk is the \
vector component of k.";

Begin["Private`"]

PolarizationSum[exp_, k_]:=Module[{uncontracted,PV,PVDD,LT,kbar},

	uncontracted=X`Utilities`Uncontract[exp, k];

	(* make internal replacements for massless pol. vectors *)
	uncontracted=uncontracted/.{LTensor[PolarizationVector[k,0],mu_]:>PV[mu]};
	uncontracted=uncontracted/.{LTensor[PolarizationVectorDag[k,0],mu_]:>PVDD[mu]};

	(* make internal replacements for massive pol. vectors *)
	uncontracted=uncontracted/.{LTensor[PolarizationVector[k,m_],mu_]:>PV[mu,m]};
	uncontracted=uncontracted/.{LTensor[PolarizationVectorDag[k,m_],mu_]:>PVDD[mu,m]};

	(* *)
	LT[x__] := LTensor[x];
	LD[x__] := LDot[x];
	kbar = Bar[k];

	ReplaceAll[uncontracted,{
		(* massless polarization spin sum *)
		PV[mu1_] * PVDD[mu2_]:> (
			-LT[MetricG,mu1,mu2]
			+ (LT[k,mu1]*LT[kbar,mu2] + LT[k,mu2]*LT[kbar,mu1]) / LD[k, kbar]
		),

		(* massive polarization spin sum *)
		PV[mu1_,m_] * PVDD[mu2_,m_]:>
			-LT[MetricG,mu1,mu2] + LT[k,mu1]*LT[k,mu2] / m^2
	}]
]

End[]
