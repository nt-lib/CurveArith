
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
    TSTAssertEQ(CAHasFunctionOfDegreeAtMost(C3, 2), false);
    TSTAssertEQ(CAHasFunctionOfDegreeAtMost(C3, 3), true);
    TSTAssertEQ(CAHasFunctionOfDegreeAtMost(C5, 3), false);
    TSTAssertEQ(CAHasFunctionOfDegreeAtMost(C5, 4), true);
end procedure;

procedure TestGonality()
    TSTAssertEQ(CAGonality(C3), 3);
    TSTAssertEQ(CAGonality(C5), 4);
end procedure;


TestHasFunctionOfDegreeAtMost();
TestGonality();