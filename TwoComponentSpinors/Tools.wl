(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

Begin["Private`"]


(* ConvertToInternal[expr]
	Convert expressions into more convinient form.
*)
ConvertToInternal[expr___]:=expr//.{
	WeylMatrixL[mtx___]:>iWML[mtx],
	WeylMatrixR[mtx___]:>iWMR[mtx],
	LTensor[WeylS, mu_] :> iLTS[mu]
};

(* ConvertToInternal[expr]
	Convert expressions back to external form.
*)
ConvertToExternal[expr___]:=expr//.{
	iWML[mtx___] :> WeylMatrixL[mtx],
	iWMR[mtx___] :> WeylMatrixR[mtx],
	iLTS[mu_] :> LTensor[WeylS, mu]
};

(* RemoveWeyl1[expr]
    Remove any unneeed Weyl1's from a WeylMatrix. Weyl1s are only necessary
    if there are not other sigma matrices.
*)
RemoveWeyl1[expr___]:=expr/.{
    (pat:WeylMatrixL|WeylMatrixR|iWML|iWMR)[left__, c_. * Weyl1,right___] :>
        c * pat[left,right],
    (pat:WeylMatrixL|WeylMatrixR|iWML|iWMR)[left___,c_. * Weyl1,right__] :>
        c * pat[left,right]
};


End[]
