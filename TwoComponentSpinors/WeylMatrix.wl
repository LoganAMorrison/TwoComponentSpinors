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

Begin["Private`"]

(* Remove any internal WeylMatrixL(R)'s from inside WeylMatrixL(R). '*)
(pat:WeylMatrixL|WeylMatrixR)[a___, (pat2:WeylMatrixL|WeylMatrixL)[b___],
    c___]:=pat[a,b,c];

End[]
