(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

(* Products of Pauli-sigma matrices *)
WeylMatrixL::usage="WeylMatrixL[\!\(\*SubscriptBox[\(mtx\), \(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)] represents the products of pauli-sigma matrices in spinor space ending with an undotted index.";

WeylMatrixR::usage="WeylMatrix[\!\(\*SubscriptBox[\(mtx\), \(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)] represents the products of pauli-sigma matrices in spinor space ending with a dotted index.";

WeylMatrixL::InvalidWeylMatrixL="WeylMatrixL[...,WeylMatrixR[...]] \
encountered. Assuming inner WeylMatrixR is incorrect and stripping it off.";

WeylMatrixR::InvalidWeylMatrixR="WeylMatrixR[...,WeylMatrixL[...]] \
encountered. Assuming inner WeylMatrixL is incorrect and stripping it off.";

Begin["Private`"]

(* Remove any internal WeylMatrixL(R)'s from inside WeylMatrixL(R). '*)

WeylMatrixL[left___, WeylMatrixR[mtx___]]:=Module[{},
    Message[WeylMatrixL::InvalidWeylMatrixL];
    WeylMatrixL[left, mtx]
];
WeylMatrixR[left___, WeylMatrixL[mtx___]]:=Module[{},
    Message[WeylMatrixR::InvalidWeylMatrixR];
    WeylMatrixR[left, mtx]
];
WeylMatrixL[left___, (WeylMatrixL|WeylMatrixR)[mtx___],right__]:=
    WeylMatrixL[left,mtx,right];
WeylMatrixR[left___, (WeylMatrixL|WeylMatrixR)[mtx___],right__]:=
    WeylMatrixL[left,mtx,right];

End[]
