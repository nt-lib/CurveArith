
P<x,y,z> := ProjectiveSpace(GF(3), 2);
f := 4*x^4 - y^4 - z^4;
// The above is an equation for the modular curve 8.96.3.e.1
// https://beta.lmfdb.org/ModularCurve/Q/8.96.3.e.1/
// This curve has gonality 3 over F_3 and gonality 4 over F_5
C3 := Curve(P,f);
P<x,y,z> := ProjectiveSpace(GF(5), 2);
f := 4*x^4 - y^4 - z^4;
C5 := Curve(P,f);

procedure TestHasFunctionOfDegreeAtMost()
    TSTAssertEQ(HasFunctionOfDegreeAtMost(C3, 2), false);
    TSTAssertEQ(HasFunctionOfDegreeAtMost(C3, 3), true);
    TSTAssertEQ(HasFunctionOfDegreeAtMost(C5, 3), false);
    TSTAssertEQ(HasFunctionOfDegreeAtMost(C5, 4), true);
end procedure;

procedure TestGonality()
    TSTAssertEQ(Gonality(C3), 3);
    TSTAssertEQ(Gonality(C5), 4);
end procedure;


TestHasFunctionOfDegreeAtMost();
TestGonality();