(* ::Package:: *)

(* Mathematica Init File *)

Needs["X`"];

BeginPackage["TwoComponentSpinors`", {"X`"}];

Print["Two Component Spinors"];

(* Definitions *)
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Weyl.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "WaveFunctions.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "WeylMatrix.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "WeylLine.wl"}]];

Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Tools.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "HermitianConjugate.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "WeylTrace.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "FermionSpinSum.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "PolarizationSum.wl"}]];

EndPackage[];
