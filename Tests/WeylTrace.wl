(* Trace of the idenitiy matrix *)
VerificationTest[WeylTrace[WeylMatrixL[]],2]
VerificationTest[WeylTrace[WeylMatrixR[]],2]


(* Trace of two sigma matrices *)
VerificationTest[
    WeylTrace[WeylMatrixL[
        LTensor[WeylS,mu],LTensor[WeylS,nu]
    ]],
    2*LTensor[MetricG,mu,nu]
]
VerificationTest[
    WeylTrace[WeylMatrixR[
        LTensor[WeylS,mu],LTensor[WeylS,nu]
    ]],
    2*LTensor[MetricG,mu,nu]
]


(* Trace of four sigma matrices *)
VerificationTest[
    WeylTrace[WeylMatrixL[
        LTensor[WeylS,mu1],
        LTensor[WeylS,mu2],
        LTensor[WeylS,mu3],
        LTensor[WeylS,mu4]
    ]],
    2*(LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]+
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4])
]
VerificationTest[
    WeylTrace[WeylMatrixR[
        LTensor[WeylS,mu1],
        LTensor[WeylS,mu2],
        LTensor[WeylS,mu3],
        LTensor[WeylS,mu4]
    ]],
    2*(LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]-
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4])
]


(* Test expression with dot products *)
VerificationTest[
    Contract[WeylTrace[WeylMatrixL[
        LDot[p1,WeylS],
        LDot[p2,WeylS],
        LDot[p3,WeylS],
        LDot[p4,WeylS]
    ]]],
    Contract[
    2*LTensor[p1,mu1]*LTensor[p2,mu2]*LTensor[p3,mu3]*LTensor[p4,mu4]*
    (LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]+
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4]
    )
    ]
]
VerificationTest[
    Contract[WeylTrace[WeylMatrixR[
        LDot[p1,WeylS],
        LDot[p2,WeylS],
        LDot[p3,WeylS],
        LDot[p4,WeylS]
    ]]],
    Contract[
    2*LTensor[p1,mu1]*LTensor[p2,mu2]*LTensor[p3,mu3]*LTensor[p4,mu4]*
    (LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]-
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4]
    )
    ]
]


(* Test expressions with sums of dot products *)
VerificationTest[
    Contract[WeylTrace[WeylMatrixL[
        LDot[p1,WeylS] + LDot[p2,WeylS],
        LDot[p3,WeylS] + LDot[p4,WeylS],
        LDot[p5,WeylS] + LDot[p6,WeylS],
        LDot[p7,WeylS] + LDot[p8,WeylS]
    ]]],
    Contract[
    2*(LTensor[p1,mu1] + LTensor[p2,mu1])*
    (LTensor[p3,mu2] + LTensor[p4,mu2])*(LTensor[p5,mu3] + LTensor[p6,mu3])*
    (LTensor[p7,mu4] + LTensor[p8,mu4])*
    (LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]+
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4]
    )
    ]
]
VerificationTest[
    Contract[WeylTrace[WeylMatrixR[
        LDot[p1,WeylS] + LDot[p2,WeylS],
        LDot[p3,WeylS] + LDot[p4,WeylS],
        LDot[p5,WeylS] + LDot[p6,WeylS],
        LDot[p7,WeylS] + LDot[p8,WeylS]
    ]]],
    Contract[
    2*(LTensor[p1,mu1] + LTensor[p2,mu1])*
    (LTensor[p3,mu2] + LTensor[p4,mu2])*(LTensor[p5,mu3] + LTensor[p6,mu3])*
    (LTensor[p7,mu4] + LTensor[p8,mu4])*
    (LTensor[MetricG,mu1,mu2]*LTensor[MetricG,mu3,mu4]-
        LTensor[MetricG,mu1,mu3]*LTensor[MetricG,mu2,mu4]+
        LTensor[MetricG,mu1,mu4]*LTensor[MetricG,mu2,mu3]-
        I*LTensor[LeviCivitaE,mu1,mu2,mu3,mu4]
    )
    ]
]
