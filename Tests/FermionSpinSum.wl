(* Test top decay to bottom + W *)
VerificationTest[
    Module[{Amp,AmpConj,MSqrdAvg},
        (* iM ~ xd_b.Sb_mu.x_t *)
        Amp=WeylLine[{1,pb,mb},{1,pt,mt},WeylMatrixL[LTensor[WeylS,mu1]]];
        AmpConj=HermitianConjugate[Amp]/.{mu1->mu2};
        Simplify[Contract[FermionSpinSum[Amp*AmpConj]]]
    ],
    Simplify[2*Contract[
        LTensor[pt,mu1]*LTensor[pb,mu2]+LTensor[pb,mu1]*LTensor[pt,mu2]-
        LTensor[MetricG,mu1,mu2]*LDot[pt,pb]-
        I*LTensor[LeviCivitaE,mu1,mu3,mu2,mu4]*LTensor[pt,mu3]*LTensor[pb,mu4]
    ]]
]


(* Test Z decay into pair of fermions *)
VerificationTest[
    Module[{Amp,AmpConj,MSqrdAvg},
        (* iM ~ af * xd_f.SB_mu.y_fb + bf * y_f.S_mu.xd_fb *)
        Amp=(
            af*WeylLine[{1,pf,mf},{-1,pfb,mf},WeylMatrixL[LTensor[WeylS,mu1]]] +
            bf*WeylLine[{-1,pf,mf},{1,pfb,mf},WeylMatrixR[LTensor[WeylS,mu1]]]
        );
        AmpConj=HermitianConjugate[Amp]/.{mu1->mu2};
        Simplify[Contract[FermionSpinSum[Amp*AmpConj]]]
    ],
    Simplify[Contract[
        af^2*WeylTrace[WeylMatrixR[
            LTensor[WeylS,mu1],LDot[pfb,WeylS],LTensor[WeylS,mu2],LDot[pf,WeylS]
        ]] +
        bf^2*WeylTrace[WeylMatrixL[
            LTensor[WeylS,mu1],LDot[pfb,WeylS],LTensor[WeylS,mu2],LDot[pf,WeylS]
        ]] -
        mf^2*af*bf*WeylTrace[WeylMatrixR[
            LTensor[WeylS,mu1],LTensor[WeylS,mu2]
        ]] -
        mf^2*af*bf*WeylTrace[WeylMatrixL[
            LTensor[WeylS,mu1],LTensor[WeylS,mu2]
        ]]
    ]]
]


(* Test muon decay into e + vmu + ve *)
VerificationTest[
    Module[{DW,Amp,AmpConj,MSqrdAvg},
        (* iM ~  *)
        DW=LDot[pmu-pvm,pmu-pvm]-MW^2;
        Amp=(-I*g/Sqrt[2])^2*(-I*LTensor[MetricG,mu1,mu2])/DW*
            WeylLine[{1,pvm,0},{1,pmu,mmu},WeylMatrixL[LTensor[WeylS,mu1]]] *
            WeylLine[{1,pe,0},{-1,pve,0},WeylMatrixL[LTensor[WeylS,mu2]]];
        AmpConj=HermitianConjugate[Amp]/.{mu1->nu1,mu2->nu2};
        Simplify[(1/2)*Contract[FermionSpinSum[Amp*AmpConj]]/.{Dim->4}]
    ],
    Simplify[
        2*g^4/(LDot[pmu-pvm,pmu-pvm]-MW^2)^2 * LDot[pe,pvm] * LDot[pmu,pve]
    ]
]
