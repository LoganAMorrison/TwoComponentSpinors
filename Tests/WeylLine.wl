(* Test that products of WeylLine's get grouped into WeylLineProduct *)
VerificationTest[
    WeylLine[{1, p1, m1}, {1, p2, m2}, WeylMatrixL[]] *
    WeylLine[{1, p3, m3}, {1, p4, m4}, WeylMatrixR[LTensor[WeylS,mu]]] *
    WeylLine[{1, p5, m5}, {1, p6, m6}, WeylMatrixL[LDot[WeylS,p]]],
    WeylLineProduct[
        WeylLine[{1, p1, m1}, {1, p2, m2}, WeylMatrixL[]],
        WeylLine[{1, p3, m3}, {1, p4, m4}, WeylMatrixR[LTensor[WeylS,mu]]],
        WeylLine[{1, p5, m5}, {1, p6, m6}, WeylMatrixL[LDot[WeylS,p]]]
    ]
]
