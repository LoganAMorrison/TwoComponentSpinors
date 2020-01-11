(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

WeylLine::usage =
  "WeylLine[{\!\(\*SubscriptBox[\(sgn\), \(1\)]\),\!\(\*SubscriptBox[\
\(p\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \
\(1\)]\)},{\!\(\*SubscriptBox[\(sgn\), \
\(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \
\(2\)]\)},WeylMatrix[\!\(\*SubscriptBox[\(mtx\), \
\(1\)]\),...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)]] represents a \
product of Weyl matrices sandwiched between two-component \
wave-function spinors \!\(\*SubscriptBox[\(\[Xi]\), \
\(1\)]\)(\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \
\(1\)]\)) and \!\(\*SubscriptBox[\(\[Xi]\), \
\(2\)]\)(\!\(\*SubscriptBox[\(p\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \
\(2\)]\)). \!\(\*SubscriptBox[\(sgn\), \(k\)]\)=1 corresponds to LH \
spinors (x or \!\(\*SuperscriptBox[\(x\), \(\[Dagger]\)]\), depending \
on structure of WeylMatrix[\!\(\*SubscriptBox[\(mtx\), \(1\)]\),...]) \
and \!\(\*SubscriptBox[\(sgn\), \(k\)]\)=-1 corresponds to RH spinors \
(y or \!\(\*SuperscriptBox[\(y\), \(\[Dagger]\)]\),depending on \
structure of WeylMatrix[...,\!\(\*SubscriptBox[\(mtx\), \(n\)]\)])";

WeylLineProduct::usage = \
"WeylLineProduct[\!\(\*SubscriptBox[\(wline\), \
\(1\)]\),\!\(\*SubscriptBox[\(wline\), \(2\)]\)] represents the \
direct product of WeylLine objects \!\(\*SubscriptBox[\(wline\), \
\(1\)]\), \!\(\*SubscriptBox[\(wline\), \(2\)]\), etc.";

Begin["`Private`"]

CurrentValue[$FrontEndSession, {InputAliases, "lx"}]= RowBox[{
        "\[LeftAngleBracket]",
        "\[ScriptX]",
        "[",
        "\[SelectionPlaceholder],",
        "\[Placeholder]",
        "]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "ly"}]= RowBox[{
        "\[LeftAngleBracket]",
        "\[ScriptY]",
        "[",
        "\[SelectionPlaceholder],",
        "\[Placeholder]",
        "]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "lxd"}]= RowBox[{
        "\[LeftAngleBracket]",
        SuperscriptBox["\[ScriptX]","\[Dagger]"],
        "[",
        "\[SelectionPlaceholder],",
        "\[Placeholder]",
        "]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "lyd"}]= RowBox[{
        "\[LeftAngleBracket]",
        SuperscriptBox["\[ScriptY]","\[Dagger]"],
        "[",
        "\[SelectionPlaceholder],",
        "\[Placeholder]",
        "]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "rx"}] = RowBox[{
    "\[ScriptX]",
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "ry"}] = RowBox[{
    "\[ScriptY]",
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "rxd"}] = RowBox[{
    SuperscriptBox["\[ScriptX]","\[Dagger]"],
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}];
CurrentValue[$FrontEndSession, {InputAliases, "ryd"}] = RowBox[{
    SuperscriptBox["\[ScriptY]","\[Dagger]"],
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}];


WeylLine /: MakeBoxes[WeylLine[{sp1:1|-1,kin1__},{sp2:1|-1,kin2__},
    (WM:WeylMatrixL|WeylMatrixR)[mtx___]],StandardForm]:= Module[{iexpr,sf1,sf2},

    iexpr=ConvertToInternal[WeylMatrixExpand[WM[mtx]]];

    If[Or[Not[FreeQ[iexpr,iWMRR]],Not[FreeQ[iexpr,iWMRL]]],
        (* If matrix has left dotted index *)
        sf1=Switch[sp1,
            1,SuperscriptBox["\[ScriptX]","\[Dagger]"],
            -1,SuperscriptBox["\[ScriptY]","\[Dagger]"]
        ],
        (* If matrix has left undotted index *)
        sf1=Switch[sp1,1,"\[ScriptX]",-1,"\[ScriptY]"];
    ];
    If[Or[Not[FreeQ[iexpr,iWMRR]],Not[FreeQ[iexpr,iWMLR]]],
        (* If matrix has right dotted index *)
        sf2=Switch[sp2,
            1,SuperscriptBox["\[ScriptX]","\[Dagger]"],
            -1,SuperscriptBox["\[ScriptY]","\[Dagger]"]
        ],
        (* If matrix has right undotted index *)
        sf2=Switch[sp2,1,"\[ScriptX]",-1,"\[ScriptY]"];
    ];

    RowBox[{"\[LeftAngleBracket]",RowBox[{sf1,"[",X`Internal`ToRowBox[{kin1}],"]"}],",",Sequence@@Riffle[Map[Function[{item},MakeBoxes[item,StandardForm],{HoldAllComplete}],Unevaluated[{mtx}]],","],",",RowBox[{sf2,"[",X`Internal`ToRowBox[{kin2}],"]"}],"\[RightAngleBracket]"}]
];


SetAttributes[WeylLineProduct,{Orderless}];

WeylLineProduct[]:=1;

WeylLineProduct[left___,WeylLineProduct[ls___],right___]:=
    WeylLineProduct[left,ls,right];

(* Change products of WeylLine's into WeylLineProduct *)
WeylLine /: Times[WeylLine[mtx1__],WeylLine[mtx2__]] :=
    WeylLineProduct[WeylLine[mtx1],WeylLine[mtx2]];

(* Change product of WeylLine with WeylLineProduct into WeylLineProduct *)
WeylLine /: Times[WeylLineProduct[lines__],WeylLine[mtx1__]] :=
    WeylLineProduct[lines,WeylLine[mtx1]];

(* Change product of WeylLineProduct's into WeylLineProduct *)
WeylLineProduct /: Times[
    WeylLineProduct[lines1__],WeylLineProduct[lines2__]] :=
    WeylLineProduct[lines1,lines2];


WeylLineProduct  /: MakeBoxes[WeylLineProduct[lines___],StandardForm]:=RowBox[List@@Riffle[Map[Function[{item},MakeBoxes[item,StandardForm],{HoldAllComplete}],Unevaluated[{lines}]],"\[CircleTimes]"]];

End[]
