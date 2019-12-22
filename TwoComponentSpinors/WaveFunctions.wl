(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: TwoComponentSpinors` *)
(* :Author: *)
(* :Date: *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 ... *)

(* Wave functions *)
SpinorX::usage="SpinorX[p, m] represents the x spinor wave-function with momentum p and mass m.";

SpinorY::usage="SpinorY[p, m] represents the y spinor wave-function with momentum p and mass m";

PolarizationVector::usage="Polarization[k, m] represents a wavefunction \ polarization vector for a vector boson with momentum k and mass m.";

Begin["Private`"]

CurrentValue[$FrontEndSession, {InputAliases, "pol"}] = RowBox[{
    "\[ScriptE]", "[",
    RowBox[{
        "\[SelectionPlaceholder]", ",", "\[Placeholder]"
        }],
        "]"
}]


End[]
