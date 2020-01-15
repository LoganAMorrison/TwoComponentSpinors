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

WeylLine::InvalidWeylMatrix="Invalid weyl-matrix encountered while \
constructing WeylLine.";

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

WeylLine /:MakeBoxes[WeylLine[{sp1 : 1 | -1, kin1__}, {sp2 : 1 | -1,kin2__}, (WM : WeylMatrixL | WeylMatrixR)[mtx___]],TraditionalForm] :=
  Module[{iexpr, sf1, sf2},
   iexpr = ConvertToInternal[WeylMatrixExpand[WM[mtx]]];
   If[Or[Not[FreeQ[iexpr, iWMRR]],
     Not[FreeQ[iexpr, iWMRL]]],(*If matrix has left dotted index*)
    sf1 = Switch[sp1, 1,
      SuperscriptBox["\[ScriptX]", "\[Dagger]"], -1,
      SuperscriptBox["\[ScriptY]",
       "\[Dagger]"]],(*If matrix has left undotted index*)
    sf1 = Switch[sp1, 1, "\[ScriptX]", -1, "\[ScriptY]"];];
   If[Or[Not[FreeQ[iexpr, iWMRR]],
     Not[FreeQ[iexpr, iWMLR]]],(*If matrix has right dotted index*)
    sf2 = Switch[sp2, 1,
      SuperscriptBox["\[ScriptX]", "\[Dagger]"], -1,
      SuperscriptBox["\[ScriptY]",
       "\[Dagger]"]],(*If matrix has right undotted index*)
    sf2 = Switch[sp2, 1, "\[ScriptX]", -1, "\[ScriptY]"];];

   RowBox[{"[", RowBox[{sf1, "(", X`Internal`ToRowBox[{kin1}], ")"}],"(",
    Sequence @@
      Riffle[
       Map[Function[{item},
         MakeBoxes[item, TraditionalForm], {HoldAllComplete}],
        Unevaluated[{mtx}]], ")("],")",
     RowBox[{sf2, "(", X`Internal`ToRowBox[{kin2}], ")"}], "]"}]
   ];

(*Parsing of Input braket notation to WeylLine object*)
MakeExpression[
    RowBox@({"\[LeftAngleBracket]",
    RowBox[{spinor2_,"[",RowBox[kin2 : {PatternSequence[_, ","] .., _}],"]"}],",",
    dMtxes : PatternSequence[_, ","] ...,
    RowBox[{spinor1_,"[",RowBox[kin1 : {PatternSequence[_, ","] .., _}],"]"}],
    "\[RightAngleBracket]"}),StandardForm]:= Module[{idx1,idx2,wm,iWM},

        idx1=If[Or[spinor2===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor2==="\[ScriptX]"],
                "1",
                RowBox[{"-","1"}]
        ];
        idx2=If[Or[spinor1===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor1==="\[ScriptX]"],
                "1",
                RowBox[{"-","1"}]
        ];
        mw=If[Or[spinor1===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor1===SuperscriptBox["\[ScriptY]","\[Dagger]"]],
              "WeylMatrixR",
              "WeylMatrixL"
        ];

        iWM=ConvertToInternal[MakeExpression[RowBox[{mw,"[",RowBox[Riffle[{dMtxes}[[1 ;; -1 ;; 2]], ","]],"]"}]]/.HoldComplete[expr___]:>expr];

        (* If iWMRR or iWMRL, need dagger on left *)
        If[Or[Head[iWM]===iWMRR,Head[iWM]===iWMRL],
            If[Not[MatchQ[spinor2,SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]]],
                Message[WeylLine::InvalidWeylMatrix];Abort[];
            ]
        ];
        (* If iWMLR or iWMLL, can't have dagger on left *)
        If[Or[Head[iWM]===iWMLR,Head[iWM]===iWMLL],
            If[Not[MatchQ[spinor2,"\[ScriptX]"|"\[ScriptY]"]],
                Message[WeylLine::InvalidWeylMatrix];Abort[];
            ]
        ];

        MakeExpression[RowBox[{"WeylLine","[",RowBox[
    	 {RowBox[{"{",RowBox[Join[{idx2,","},Riffle[kin2[[1 ;; -1 ;; 2]], ","]]],"}"}],",",
    	  RowBox[{"{",RowBox[Join[{idx1,","},Riffle[kin1[[1 ;; -1 ;; 2]], ","]]],"}"}],",",
    	  RowBox[{mw,"[",RowBox[Riffle[{dMtxes}[[1 ;; -1 ;; 2]], ","]],"]"}]
    	 }],"]"}],StandardForm
    	]
]/;And[MatchQ[spinor2,"\[ScriptX]"|"\[ScriptY]"|SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]],
       MatchQ[spinor1,"\[ScriptX]"|"\[ScriptY]"|SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]]];

MakeExpression[
    RowBox@({"\[LeftAngleBracket]",RowBox[{RowBox[{spinor2_,"[",RowBox[kin2 : {PatternSequence[_, ","] .., _}],"]"}],",",
    dMtxes : PatternSequence[_, ","] ...,
    RowBox[{spinor1_,"[",RowBox[kin1 : {PatternSequence[_, ","] .., _}],"]"}]}],"\[RightAngleBracket]"}),StandardForm]:=Module[{idx1,idx2,wm,iWM},

        idx2=If[Or[spinor2===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor2==="\[ScriptX]"],
                "1",
                RowBox[{"-","1"}]
        ];
        idx1=If[Or[spinor1===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor1==="\[ScriptX]"],
                "1",
                RowBox[{"-","1"}]
        ];
        mw=If[Or[spinor1===SuperscriptBox["\[ScriptX]","\[Dagger]"], spinor1===SuperscriptBox["\[ScriptY]","\[Dagger]"]],
              "WeylMatrixR",
              "WeylMatrixL"
        ];

        iWM=ConvertToInternal[MakeExpression[RowBox[{mw,"[",RowBox[Riffle[{dMtxes}[[1 ;; -1 ;; 2]], ","]],"]"}]]/.HoldComplete[expr___]:>expr];

        (* If iWMRR or iWMRL, need dagger on left *)
        If[Or[Head[iWM]===iWMRR,Head[iWM]===iWMRL],
            If[Not[MatchQ[spinor2,SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]]],
                Message[WeylLine::InvalidWeylMatrix];Abort[];
            ]
        ];
        (* If iWMLR or iWMLL, can't have dagger on left *)
        If[Or[Head[iWM]===iWMLR,Head[iWM]===iWMLL],
            If[Not[MatchQ[spinor2,"\[ScriptX]"|"\[ScriptY]"]],
                Message[WeylLine::InvalidWeylMatrix];Abort[];
            ]
        ];

        MakeExpression[RowBox[{"WeylLine","[",RowBox[
    	 {RowBox[{"{",RowBox[Join[{idx2,","},Riffle[kin2[[1 ;; -1 ;; 2]], ","]]],"}"}],",",
    	  RowBox[{"{",RowBox[Join[{idx1,","},Riffle[kin1[[1 ;; -1 ;; 2]], ","]]],"}"}],",",
    	  RowBox[{mw,"[",RowBox[Riffle[{dMtxes}[[1 ;; -1 ;; 2]], ","]],"]"}]
    	 }],"]"}],StandardForm
    	]
]/;And[MatchQ[spinor2,"\[ScriptX]"|"\[ScriptY]"|SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]],
       MatchQ[spinor1,"\[ScriptX]"|"\[ScriptY]"|SuperscriptBox["\[ScriptX]","\[Dagger]"]|SuperscriptBox["\[ScriptY]","\[Dagger]"]]];


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
WeylLineProduct  /: MakeBoxes[WeylLineProduct[lines___],TraditionalForm]:=RowBox[List@@Riffle[Map[Function[{item},MakeBoxes[item,TraditionalForm],{HoldAllComplete}],Unevaluated[{lines}]],"\[CircleTimes]"]];

End[]
