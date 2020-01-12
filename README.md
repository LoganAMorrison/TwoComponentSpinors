# TwoComponentSpinors
Mathematica extension of Package-X for two-component spinors.

## Installation
Clone or download this repository into your Mathematica applications directory.
Then, running

```Mathematica
<<TwoComponentSpinors`
```

should load Package-X and the `TwoComponentSpinors` extension.

## Usage

### Basics

There are two types of matrices in `TwoComponentSpinors`: ones that end with
an *undotted* spinor index and ones that end with a *dotted* spinor index. Matrices that end with and *undotted* spinor index are called `WeylMatrixL` and those that end with a *dotted* spinor index are called `WeylMatrixR`. For example:

```Mathematica
(* represents sigma_mu * sigmabar_mu *)
WeylMatrixL[LTensor[WeylS,mu],LTensor[WeylS,nu]]

(* represents sigmabar_mu * sigma_mu *)
WeylMatrixR[LTensor[WeylS,mu],LTensor[WeylS,nu]]
```

Here, `WeylS` represents the Pauli sigma 4-vector: WeylS_mu=(1,sigma_1,sigma_2,sigma_3).

### Traces of sigma matrices


```Mathematica
(* Two sigma matrices *)

WeylTrace[WeylMatrixL[LTensor[WeylS, mu], LTensor[WeylS, nu]]]
(* Result: 2*LTensor[MetricG, mu, nu] *)

WeylTrace[WeylMatrixR[LTensor[WeylS, mu], LTensor[WeylS, nu]]]
(* Result: 2*LTensor[MetricG, mu, nu] *)


(* Four sigma matrices *)
WeylTrace[WeylMatrixL[
    LTensor[WeylS, mu1], LTensor[WeylS, mu2],
    LTensor[WeylS, mu3], LTensor[WeylS, mu4]
]]
(* Result:
    2*(LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu3] -
       LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu4] +
       LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu4] +
       I*LTensor[LeviCivitaE, mu1, mu2, mu3, mu4])
*)

WeylTrace[WeylMatrixR[
    LTensor[WeylS, mu1], LTensor[WeylS, mu2],
    LTensor[WeylS, mu3], LTensor[WeylS, mu4]
]]
(* Result:
    2*(LTensor[MetricG, mu1, mu4]*LTensor[MetricG, mu2, mu3] -
       LTensor[MetricG, mu1, mu3]*LTensor[MetricG, mu2, mu4] +
       LTensor[MetricG, mu1, mu2]*LTensor[MetricG, mu3, mu4] +
       -I*LTensor[LeviCivitaE, mu1, mu2, mu3, mu4])
*)
```

### Spinor Lines and  Hermitian Conjugates

Construct a spinor line with a massive polarization vector:

```Mathematica
exp = WeylLine[{1,p1,m1},{1,p2,m2}. WeylMatrixR[LTensor[WeylS, mu]]]
       * LTensor[PolarizationVectorDag[k, mz], mu];
```
