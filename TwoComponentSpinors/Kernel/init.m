(* ::Package:: *)

(* Mathematica Init File *)

Needs["X`"];

BeginPackage["TwoComponentSpinors`"];

Print["Two Component Spinors"];

Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Pauli.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "WaveFunctions.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "SigmaMatrix.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "SpinorLine.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "HermitianConjugate.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Uncontract.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "SpinorTrace.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "FermionSpinSum.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "PolarizationSum.m"}]];

EndPackage[];
