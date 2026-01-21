procedure TestClassGroupCrv()
    P<x,y,z> := ProjectiveSpace(GF(3), 2);
    f := 4*x^4 - y^4 - z^4;
    C := Curve(P,f);
    CG := ClassGroup(C);
    TSTAssertEQ(AbelianInvariants(CG), [4, 4, 4, 0]);
end procedure;


procedure TestClassGroupFldFunG()
    QQ := Rationals();
    P<x,y> := PolynomialRing(QQ, 2);
    f := 4*x^4 - y^4 - 1;
    /* the above is an equation for the modular curve 8.96.3.e.1
    https://beta.lmfdb.org/ModularCurve/Q/8.96.3.e.1/
    the Jacobian of this curve is isogenous to the E^3 where E is an 
    elliptic curve in the isogeny class 64.a
    https://beta.lmfdb.org/EllipticCurve/Q/64/a/
    the traces of Frobenius of E are given by:
    */
    traces := [
        [3,  0],
        [5,  2],
        [7,  0],
        [11, 0],
        [13, -6],
        [17,  2],
        [19,  0],
        [23,  0],
        [29, 10]
    ];
    for pair in traces do
        p := pair[1];
        print p;
        trace := pair[2];
        Kp := FunctionField(ChangeRing(f, GF(p)));
        time CG := ClassGroup(Kp);
        TSTAssertEQ(TorsionFreeRank(CG), 1);
        TSTAssertEQ(Order(TorsionSubgroup(CG)), (1 + p - trace)^3);
    end for;
end procedure;

/*
TestClassGroupCrv();
the above currently failes with the following error:
Loading "src/test_classgroup.m"

TestClassGroup(
)
ClassGroup(
    C: C
)
In file "/opt/magma2.28-3/package/Geometry/Crv/class_group.m", line 13, column 
16:
>>     Cl, _, phi := ClassGroup(FF);
                  ^
Runtime error in :=: Expected to assign 3 value(s) but only computed 1 value(s)
*/

TestClassGroupFldFunG();