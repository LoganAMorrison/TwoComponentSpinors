(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

(*WeylLineProductExpand::usage="";*)

Begin["`Private`"]


(* ConvertToInternal[expr]
	Convert expressions into more convinient form.
*)
ConvertToInternal[expr___]:=expr//.{
    WeylMatrixR[c_. Weyl1] :> c * iWMRR[],
    WeylMatrixL[c_. Weyl1] :> c * iWMLL[],
	LTensor[WeylS, mu_] :> iLTS[mu],
    WeylMatrixR[mtx___] :> If[OddQ[Count[List[mtx],pat:(c_. iLTS[_] + b_.)|(c_. LTensor[WeylS, _] + b_.)|(c_. LDot[WeylS, _] + b_.)]],iWMLR[mtx],iWMRR[mtx]],
	WeylMatrixL[mtx___] :> If[OddQ[Count[List[mtx],pat:(c_. iLTS[_] + b_.)|(c_. LTensor[WeylS, _] + b_.)|(c_. LDot[WeylS, _] + b_.)]],iWMRL[mtx],iWMLL[mtx]]
};

(* ConvertToInternal[expr]
	Convert expressions back to external form.
*)
ConvertToExternal[expr___]:=expr//.{
	iLTS[mu_] :> LTensor[WeylS, mu],
    (iWMRR|iWMLR)[mtx___] :> WeylMatrixR[mtx],
    (iWMLL|iWMRL)[mtx___] :> WeylMatrixL[mtx]
};

(* RemoveWeyl1[expr]
    Remove any Weyl1's from a WeylMatrix.
*)
RemoveWeyl1[expr___]:=expr//.{
    (pat:WeylMatrixL|WeylMatrixR|iWMLL|iWMRR|iWMLR|iWMRL)[left___, c_. * Weyl1,right___] :> c * pat[left,right]
};

(* WeylMatrixExpand[WeylMatrixL[mtx]] or WeylMatrixExpand[WeylMatrixR[mtx]]
    Use linearity of matrix multiplication to expand out Weyl matrices.
*)
WeylMatrixExpand[(pat:WeylMatrixL|WeylMatrixR)[mtx___]]:=Module[{iexpr},
    iexpr=ExpandAll[pat[mtx]];
    iexpr=Distribute[iexpr,Plus,pat];
    iexpr=X`Utilities`Uncontract[iexpr,WeylS];
    (* Remove Constants *)
    iexpr=iexpr//.{
        pat[left___, c_ LTensor[WeylS, mu_], right___] :> c * pat[left, LTensor[WeylS, mu], right],
        pat[left___, c_. Weyl1, right___] :> c * pat[left, right]
    }
];


(* WeylLineExpand[WeylLine[mtx]]
    Use linearity of matrix multiplication to expand out Weyl matrices.
*)
WeylLineExpand[WeylLine[kin1_List, kin2_List, (pat:WeylMatrixL|WeylMatrixR)[mtx___]]]:=Module[{iexpr},
    iexpr=WeylLine[kin1, kin2, WeylMatrixExpand[pat[mtx]]];
    iexpr=Distribute[iexpr, Plus, WeylLine];
    iexpr=iexpr//.{
        WeylLine[kin1, kin2, c_ pat[imtx___]] :> c * WeylLine[kin1, kin2, pat[imtx]]
    }
];

WeylLineProductExpand[WeylLineProduct[lines___]]:=Module[{iexpr},
    (* Expand out inner WeylLines *)
    iexpr=ExpandAll[WeylLineProduct[lines]/.{
        WeylLine[kin1_List, kin2_List, (pat:WeylMatrixL|WeylMatrixR)[mtx___]]:> WeylLineExpand[WeylLine[kin1,kin2,pat[mtx]]]
    },WeylLine[___]];
    iexpr=Distribute[iexpr, Plus, WeylLineProduct];
    iexpr=iexpr//.{
        WeylLineProduct[left___, c_ WeylLine[line___], right___] :> c * WeylLineProduct[left, WeylLine[line], right]
    };
    Collect[iexpr, WeylLineProduct[___]]
];


End[]
