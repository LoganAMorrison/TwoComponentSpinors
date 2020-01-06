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
WeylS::usage="WeylS represents the Pauli matrix in two-component spinor space";
Weyl1::usage="Weyl1 represents the unit matrix in two-component spinor space";

Begin["Private`"]

(*
WeylLine /: MakeBoxes[Weyl1,StandardForm]:=

    RowBox[{"\[LeftAngleBracket]",RowBox[{sf1,"[",X`Internal`ToRowBox[{kin1}],"]"}],",",Sequence@@Riffle[Map[
        Function[{item},
            MakeBoxes[item,StandardForm],
            {HoldAllComplete}
        ],Unevaluated[{mtx}]],","],",",RowBox[{sf2,"[",X`Internal`ToRowBox[{kin2}],"]"}],"\[RightAngleBracket]"}]
];
*)


End[]

EndPackage[]
