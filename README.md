# TwoComponentSpinors
Mathematica extension of PackageX for two-component spinors

## Usage

### Traces of sigma matrices

```Mathematica
(* Two sigma matrices *)
SpinorTrace[SigmaMatrix[LTensor[PauliSBar, mu], LTensor[PauliS, nu]]]
(* Result: 2*LTensor[MetricG, mu, nu] *)

(* Four sigma matrices *)
SpinorTrace[SigmaMatrix[
    LTensor[PauliSBar, mu1],
    LTensor[PauliS, mu2],
    LTensor[PauliSBar, mu3],
    LTensor[PauliS, mu4]
]]
(* Result:
    2*(LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu3]
        - LTensor[MetricG, mu1, mu3]* LTensor[MetricG, mu2, mu4]
        + LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu4]
        - I*LTensor[LeviCivitaE, mu1, mu2, mu3, mu4])
*)

(* Six sigma matrices *)
Contract[SpinorTrace[SigmaMatrix[
    LTensor[PauliS, mu1]
    LTensor[PauliSBar, mu2],
    LTensor[PauliS, mu3],
    LTensor[PauliSBar, mu4],
    LTensor[PauliS, mu5],
    LTensor[PauliSBar, mu6]
]]]
(* Result:
    -6*LTensor[MetricG, mu1, mu6]*LTensor[MetricG, mu2, mu5]*LTensor[MetricG, mu3, mu4] +
 2*Dim*LTensor[MetricG, mu1, mu6]*LTensor[MetricG, mu2, mu5]*LTensor[MetricG, mu3, mu4] +
 6*LTensor[MetricG, mu1, mu5]*LTensor[MetricG, mu2, mu6]*LTensor[MetricG, mu3, mu4] -
 2*Dim*LTensor[MetricG, mu1, mu5]*LTensor[MetricG, mu2, mu6]*LTensor[MetricG, mu3, mu4] +
 6*LTensor[MetricG, mu1, mu6]*LTensor[MetricG, mu2, mu4]*LTensor[MetricG, mu3, mu5] -
 2*Dim*LTensor[MetricG, mu1, mu6]*LTensor[MetricG, mu2, mu4]*LTensor[MetricG, mu3, mu5] -
 6*LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu6]*LTensor[MetricG, mu3, mu5] +
 2*Dim*LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu6]*LTensor[MetricG, mu3, mu5] -
 6*LTensor[MetricG, mu1, mu5]*LTensor[MetricG, mu2, mu4]*LTensor[MetricG, mu3, mu6] +
 2*Dim*LTensor[MetricG, mu1, mu5]*LTensor[MetricG, mu2, mu4]*LTensor[MetricG, mu3, mu6] +
 6*LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu5]*LTensor[MetricG, mu3, mu6] -
 2*Dim*LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu5]*LTensor[MetricG, mu3, mu6] +
 2*LTensor[MetricG, mu1, mu6]*LTensor[MetricG, mu2, mu3]*LTensor[MetricG, mu4, mu5] -
 2*LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu6]*LTensor[MetricG, mu4, mu5] +
 2*LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu6]*LTensor[MetricG, mu4, mu5] -
 2*LTensor[MetricG, mu1, mu5]*LTensor[MetricG, mu2, mu3]*LTensor[MetricG, mu4, mu6] +
 2*LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu5]*LTensor[MetricG, mu4, mu6] -
 2*LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu5]*LTensor[MetricG, mu4, mu6] +
 2*LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu3]*LTensor[MetricG, mu5, mu6] -
 2*LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu4]*LTensor[MetricG, mu5, mu6] +
 2*LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu4]*LTensor[MetricG, mu5, mu6] +
 (2*I)*LTensor[MetricG, mu5, mu6]*LTensor[LeviCivitaE, mu1, mu2, mu3, mu4] -
 (2*I)*LTensor[MetricG, mu4, mu6]*LTensor[LeviCivitaE, mu1, mu2, mu3, mu5] +
 (2*I)*LTensor[MetricG, mu4, mu5]*LTensor[LeviCivitaE, mu1, mu2, mu3, mu6] +
 (2*I)*LTensor[MetricG, mu2, mu3]*LTensor[LeviCivitaE, mu1, mu4, mu5, mu6] -
 (2*I)*LTensor[MetricG, mu1, mu3]*LTensor[LeviCivitaE, mu2, mu4, mu5, mu6] +
 (2*I)*LTensor[MetricG, mu1, mu2]*LTensor[LeviCivitaE, mu3, mu4, mu5, mu6]
*)
```

### Spinor Lines and  Hermitian Conjugates

Construct a spinor line with a massive polarization vector:

```Mathematica
exp = SpinorLine[
    SpinorXDag[p2, m2],
    SigmaMatrix[LTensor[PauliSBar, mu]],
    SpinorX[p1, m1]
    ]* LTensor[PolarizationVectorDag[k, mz], mu];
```

Compute the hermitian conjugate:

```Mathematica
expdag = HermitianConjugate[exp] /. {mu -> nu}
(* Result:
    LTensor[PolarizationVector[k, mz], nu] * SpinorLine[
        SpinorXDag[p1, m1],
        SigmaMatrix[LTensor[PauliSBar, nu]],
        SpinorX[p2, m2]
    ]
*)
```

### Computing Fermion Spin Sums and Polarization Sums

Compute the fermion spin sum and polarization sum of the above example:

```Mathematica
Contract[PolarizationSum[FermionSpinSum[exp*expdag]]]
(* Result:
    (4*LDot[k, p1]*LDot[k, p2])/mz^2
    - 4*LDot[p1, p2] + 2*Dim*LDot[p1, p2]
    - (2*LDot[k, k]*LDot[p1, p2])/mz^2
*)
```


```
