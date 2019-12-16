(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

BeginPackage["TwoComponentSpinors`"]

FermionSpinSum::usage="";

Begin["Private`"]

FermionSpinSum[exp_]:=Module[{ret},
	exp/.{
		SpinorLine[
            SpinorX[p1_,m1_],
            SigmaMatrix[sig1__],
            SpinorX[p2_,m2_]
        ] *
        SpinorLine[
            SpinorXDag[p2_,m2_],
            SigmaMatrix[sig2__],
            SpinorXDag[p1_,m1_]
        ]:>SpinorTrace[
            SigmaMatrix[sig1, X`LDot[p2,PauliS], sig2, X`LDot[p1,PauliSBar]]],
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


End[]

EndPackage[]
