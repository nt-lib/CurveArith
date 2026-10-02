
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

// A hyperelliptic curve over F_3 without rational places; x has degree 2 and its poles form
// a single place of degree 2, so it cannot be padded up to a divisor of degree 3.
K<x> := FunctionField(GF(3));
R<y> := PolynomialRing(K);
F := FunctionField(y^2 - (2*x^8 + x^2 + 2));

procedure TestHasFunctionOfDegreeAtMostNoRationalPlaces()
    TSTAssertEQ(#Places(F, 1), 0);
    for Al in ["LinAlg", "Hess"] do
        TSTAssertEQ(CAHasFunctionOfDegreeAtMost(F, 1 : Al := Al), false);
        TSTAssertEQ(CAHasFunctionOfDegreeAtMost(F, 2 : Al := Al), true);
        TSTAssertEQ(CAHasFunctionOfDegreeAtMost(F, 3 : Al := Al), true);
        TSTAssertEQ(CAGonality(F : Al := Al), 2);
    end for;
end procedure;


TestHasFunctionOfDegreeAtMost();
TestGonality();
TestHasFunctionOfDegreeAtMostNoRationalPlaces();