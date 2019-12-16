(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

SpinorTrace::usage="";

Begin["Private`"]


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


End[]

EndPackage[]
