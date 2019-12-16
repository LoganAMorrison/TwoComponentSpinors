(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

(* Pauli matrices *)
PauliS::usage="PauliS represents the Pauli matrix in two-component spinor space";
PauliSBar::usage="PauliSBar represents the barred Pauli matrix in two-component spinor space";
Pauli1::usage="Pauli1 represents the unit matrix in two-component spinor space";

(* Wave functions *)
SpinorX::usage="SpinorX[p, m] represents the x spinor wave-function with momentum p and mass m";
SpinorXDag::usage="SpinorXDag[p, m] represents the \!\(\*SuperscriptBox[\(x\), \(\[Dagger]\)]\) spinor wave-function with momentum p and mass m";
SpinorY::usage="SpinorY[p, m] represents the y spinor wave-function with momentum p and mass m";
SpinorYDag::usage="SpinorYDag[p, m] represents the \!\(\*SuperscriptBox[\(y\), \(\[Dagger]\)]\) spinor wave-function with momentum p and mass m";
PolarizationVector::usage="Polarization[k, m]";
PolarizationVectorDag::usage="PolarizationVectorDag[k,m]";

(* Products of Pauli-sigma matrices *)
SigmaMatrix::usage="SigmaMatrix[\!\(\*SubscriptBox[\(mtx\), \(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)] represents the products of pauli-sigma matrices in spinor space.";

(* Spinor Lines *)
SpinorLine::usage="SpinorLine[wf1,SigmaMatrix[...],wf2] represents are product of \
pauli matrices sandwiched between wave-functions wf1 and wf2.";

HermitianConjugate::usage="HermitianConjugate[x] returns the hermitian conjugate of x.";

SpinorTrace::usage="";
SpinSum::usage="";
PolarizationSum::usage="";



Begin["Private`"]

HermitianConjugate[PolarizationVector[p_,m_]]:=PolarizationVectorDag[p,m];
HermitianConjugate[PolarizationVectorDag[p_,m_]]:=PolarizationVector[p,m];

HermitianConjugate[SpinorX[p_,m_]]:=SpinorXDag[p,m];
HermitianConjugate[SpinorXDag[p_,m_]]:=SpinorX[p,m];

HermitianConjugate[SpinorY[p_,m_]]:=SpinorYDag[p,m];
HermitianConjugate[SpinorYDag[p_,m_]]:=SpinorY[p,m];

HermitianConjugate[SigmaMatrix[a__]]:=SigmaMatrix[Reverse[List[a]]]/.{List[b__]:>b};

HermitianConjugate[SpinorLine[a_,b_,c_]]:= SpinorLine[HermitianConjugate[c],HermitianConjugate[b],HermitianConjugate[a]];

HermitianConjugate[exp_]:= Module[{expanded},
	ReplaceAll[Expand[exp],{
		PolarizationVector[p_,m_]:>HermitianConjugate[PolarizationVector[p,m]],
		PolarizationVectorDag[p_,m_]:>HermitianConjugate[PolarizationVectorDag[p,m]],
		SpinorX[p_,m_]:>HermitianConjugate[SpinorX[p,m]],
		SpinorXDag[p_,m_]:>HermitianConjugate[SpinorXDag[p,m]],
		SpinorY[p_,m_]:>HermitianConjugate[SpinorY[p,m]],
		SpinorYDag[p_,m_]:>HermitianConjugate[SpinorYDag[p,m]],
		SigmaMatrix[a__]:>HermitianConjugate[SigmaMatrix[a]],
		SpinorLine[a_,b_,c_]:>HermitianConjugate[SpinorLine[a,b,c]]
	}]
];

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

SpinorTrace::InvalidSigmaMatrix="Invalid SigmaMatrix input in SpinorTrace[]";
SpinorTrace[exp_]:=Module[{args,evenargs,oddargs,sigmtx,coeff,mu},
	(* uncontract all arguments *)
	{coeff,sigmtx}=FixedPoint[UncontractSigmaMatrix,exp]/.{c_. SigmaMatrix[a__]:>{c,SigmaMatrix[a]}};

	(* grab all the arguments of SigmaMatrix *)
	args=sigmtx/.SigmaMatrix[x__]:>List[x];

	If[OddQ[Length[args]],Throw["InvalidSigmaMatrix",SpinorTrace::InvalidSigmaMatrix]];

	(* perform checks on arguments *)
	evenargs=Table[args[[2*i]]/.X`LTensor[x_,y_]:>x,{i,Floor[Length[args]/2]}];
	oddargs=Table[args[[2*i-1]]/.X`LTensor[x_,y_]:>x,{i,Ceiling[Length[args]/2]}];
	If[Or[Length[DeleteDuplicates[evenargs]]!=1,Length[DeleteDuplicates[oddargs]]!=1],Throw["InvalidSigmaMatrix",SpinorTrace::InvalidSigmaMatrix]];

	(* otherwise, recurrsively evaluate *)
	coeff * ReplaceRepeated[sigmtx,{
	(* base case of two sigma matrices *)
	SigmaMatrix[ X`LTensor[PauliS, mu1_], X`LTensor[PauliSBar, mu2_] ] :> 2 * X`LTensor[X`MetricG, mu1, mu2],
	SigmaMatrix[ X`LTensor[PauliSBar, mu1_], X`LTensor[PauliS, mu2_] ] :> 2 * X`LTensor[X`MetricG, mu1, mu2],
	(* base case of four sigma matrices *)
	SigmaMatrix[ X`LTensor[PauliS, mu1_], X`LTensor[PauliSBar, mu2_], X`LTensor[PauliS, mu3_], X`LTensor[PauliSBar, mu4_] ] :> 2 * (
		X`LTensor[X`MetricG, mu1, mu2]*X`LTensor[X`MetricG, mu3, mu4] -
		X`LTensor[X`MetricG, mu1, mu3]*X`LTensor[X`MetricG, mu2, mu4] +
		X`LTensor[X`MetricG, mu1, mu4]*X`LTensor[X`MetricG, mu2, mu3] +
		I * X`LTensor[X`LeviCivitaE, mu1, mu2, mu3, mu4]
	),
	SigmaMatrix[ X`LTensor[PauliSBar, mu1_], X`LTensor[PauliS, mu2_], X`LTensor[PauliSBar, mu3_], X`LTensor[PauliS, mu4_] ] :> 2 * (
		X`LTensor[X`MetricG, mu1, mu2]*X`LTensor[X`MetricG, mu3, mu4] -
		X`LTensor[X`MetricG, mu1, mu3]*X`LTensor[X`MetricG, mu2, mu4] +
		X`LTensor[X`MetricG, mu1, mu4]*X`LTensor[X`MetricG, mu2, mu3] -
		I * X`LTensor[X`LeviCivitaE, mu1, mu2, mu3, mu4]
	),
	SigmaMatrix[ X`LTensor[PauliS, mu1_], X`LTensor[PauliSBar, mu2_], X`LTensor[PauliS, mu3_], a__ ] :> (
		X`LTensor[X`MetricG, mu1, mu2] * SigmaMatrix[ X`LTensor[PauliS, mu3], a ] -
		X`LTensor[X`MetricG, mu1, mu3] * SigmaMatrix[ X`LTensor[PauliS, mu2], a ] +
		X`LTensor[X`MetricG, mu2, mu3] * SigmaMatrix[ X`LTensor[PauliS, mu1], a ] +
		I * X`LTensor[X`LeviCivitaE, mu1, mu2, mu3, mu] * SigmaMatrix[ X`LTensor[PauliS, mu], a ]
	),
	SigmaMatrix[ X`LTensor[PauliSBar, mu1_], X`LTensor[PauliS, mu2_], X`LTensor[PauliSBar, mu3_], a__ ] :> (
		X`LTensor[X`MetricG, mu1, mu2] * SigmaMatrix[ X`LTensor[PauliSBar, mu3], a ] -
		X`LTensor[X`MetricG, mu1, mu3] * SigmaMatrix[ X`LTensor[PauliSBar, mu2], a ] +
		X`LTensor[X`MetricG, mu2, mu3] * SigmaMatrix[ X`LTensor[PauliSBar, mu1], a ] -
		I * X`LTensor[X`LeviCivitaE, mu1, mu2, mu3, mu] * SigmaMatrix[ X`LTensor[PauliSBar, mu], a ]
	)
}]

];


SpinSum[exp_]:=Module[{ret},
	exp/.{
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliS],sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*SpinorTrace[SigmaMatrix[sig1,X`LDot[p2,PauliSBar],sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorX[p2_,m2_]]*SpinorLine[SpinorY[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorXDag[p2_,m2_]]*SpinorLine[SpinorYDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>-m2*m1*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>m2*m1*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>m2*m1*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorY[p2_,m2_]]*SpinorLine[SpinorX[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>-m2*m1*SpinorTrace[SigmaMatrix[sig1,sig2]],

		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliSBar]]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>m2*SpinorTrace[SigmaMatrix[sig1,sig2,X`LDot[p1,PauliS]]],
		SpinorLine[SpinorY[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorX[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorYDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorXDag[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorX[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorY[p1_,m1_]]:>-m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]],
		SpinorLine[SpinorXDag[p1_,m1_],SigmaMatrix[sig1__],SpinorYDag[p2_,m2_]]*SpinorLine[SpinorXDag[p2_,m2_],SigmaMatrix[sig2__],SpinorYDag[p1_,m1_]]:>m1*m2*SpinorTrace[SigmaMatrix[sig1,sig2]]
	}
];


PolarizationSum[exp_]:=Module[{uncontracted},
	uncontracted=FixedPoint[UncontractPolarizationVector,exp];
	ReplaceAll[uncontracted,{
		X`LTensor[PolarizationVector[k_,0],mu1_]*X`LTensor[PolarizationVectorDag[k_,0],mu2_]:>-X`LTensor[X`MetricG,mu1,mu2],
		X`LTensor[PolarizationVector[k_,m_],mu1_]*X`LTensor[PolarizationVectorDag[k_,m_],mu2_]:>-X`LTensor[X`MetricG,mu1,mu2]+X`LTensor[k,mu1]X`LTensor[k,mu2]/m^2
	}]
]



End[]

EndPackage[]
