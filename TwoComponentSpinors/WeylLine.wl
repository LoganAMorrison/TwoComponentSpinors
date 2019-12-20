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

Begin["Private`"]

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
}]
CurrentValue[$FrontEndSession, {InputAliases, "ry"}] = RowBox[{
    "\[ScriptY]",
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}]
CurrentValue[$FrontEndSession, {InputAliases, "rxd"}] = RowBox[{
    SuperscriptBox["\[ScriptX]","\[Dagger]"],
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}]
CurrentValue[$FrontEndSession, {InputAliases, "ryd"}] = RowBox[{
    SuperscriptBox["\[ScriptY]","\[Dagger]"],
    "[",
    "\[SelectionPlaceholder]",
    ",",
    "\[Placeholder]",
    "]",
    "\[RightAngleBracket]"
}]




(*
WeylLine /: MakeBoxes[
    WeylLine[{sp2:1|-1,kin2__},{sp1:1|1,kin1__},
    WeylMatrix[mtx__]],StandardForm]:=
	  With[{sf1=Switch[sp1,1,"\[ScriptU]",-1,"\[ScriptV]"],
			sf2=Switch[sp2,1,"\[ScriptU]",-1,"\[ScriptV]"]},
		RowBox[{"\[LeftAngleBracket]",RowBox[{sf2,"[",X`Internal`ToRowBox[{kin2}],"]"}],",",Sequence@@Riffle[Map[Function[{item},MakeBoxes[item,StandardForm],{HoldAllComplete}],Unevaluated[{mtx}]],","],",",RowBox[{sf1,"[",X`Internal`ToRowBox[{kin1}],"]"}],"\[RightAngleBracket]"}]
	  ];

FermionLine /: MakeBoxes[FermionLine[{sp2:1|-1,kin2__},{sp1:1|-1,kin1__},DiracMatrix[]],StandardForm]:=
	  With[{sf1=Switch[sp1,1,"\[ScriptU]",-1,"\[ScriptV]"],
			sf2=Switch[sp2,1,"\[ScriptU]",-1,"\[ScriptV]"]},
		RowBox[{"\[LeftAngleBracket]",RowBox[{sf2,"[",X`Internal`ToRowBox[{kin2}],"]"}],",",RowBox[{sf1,"[",X`Internal`ToRowBox[{kin1}],"]"}],"\[RightAngleBracket]"}]
	  ];
*)


End[]
