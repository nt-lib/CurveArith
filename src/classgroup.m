declare verbose ClassGroup, 1;

declare attributes FldFunG: classgroup_data;

classgroup_data_format := recformat<
    factor_basis,           // Factor basis
    factor_basis_count,     // Number of places in the factor basis of each degree
    factor_basis_minima,    // The indices of the places of the factor basis grouped by the underlying place in k(t)
    infinite_places,        // The indices of the infinite places of the factor basis
    relations,              // Matrix of relations
    empty_columns,          // Set of empty columns of the relation matrix
    rank,                   // Rank of relation matrix
    h,                      // Class number
    h_bounds,               // Bounds on the class number
    elementary_divisors,    // Structure of the class group
    is_complete,            // Whether the stored matrix of relations describes the full relation lattice
    log                     // Logging data
>;

log_format := recformat<
    relation_time,          // Total time taken for relation finding
    rank_time,              // Total time taken for rank computations
    h_time,                 // Total time taken for elementary divisor computations
    functions_tried,        // Total number of functions tried
    relations_found,        // Total number of relations found
    latest_relation_time,   // Time taken for most recent relation finding iteration
    latest_rank_time,       // Time taken for most recent rank computation
    latest_h_time,          // Time taken for most recent elementary divisor computation
    latest_functions_tried, // Number of functions tried in most recent relation finding iteration
    latest_relations_found  // Number of relations found in most recent relation finding iteration
>;

function NumberOfSmoothDivisors(n, m, P) 
/*
Given two integers n and m and a sequence P containing the number of places of degree <= min(n, m) on a function field F,
returns the number of effective divisors on F of degree n consisting only of places of degree <= m.
The formula used is due to Hess.
*/
    if n lt 0 or m lt 0 then
        return 0;
    end if;
    if n eq 0 then
        return 1;
    end if;
    assert #P ge Minimum(n, m);

    N := [&+[d * P[d] : d in Divisors(i) | d le m] : i in [1..n]];
    psi := [1]; // First entry is psi(0, m)
    for i in [1..n] do
        Append(~psi, &+[N[i - j] * psi[j + 1] : j in [0..i-1]] div i);
    end for;
    
    return psi[n + 1];
end function;

function LeadingCoefficients(fs, p)
/*
Given a sequence of functions fs in a function field, and a place p on the same function field,
returns nth coefficient of the Laurent series expansion of the fs at p, where -n is the
minimum valuation of all fs at p.
*/
    if #fs eq 0 then
        return [];
    end if;

    vD := Minimum([Valuation(f, p) : f in fs]);
    if vD eq 0 then
        return [Evaluate(f, p) : f in fs];
    end if;
    
    coeff_field := ResidueClassField(p);
    evaluations := [];
    for f in fs do
        if Valuation(f, p) gt vD then
            Append(~evaluations, coeff_field!0);
        else
            if assigned uniformizer then
                Append(~evaluations, Evaluate(f/uniformizer, p));
            else
                uniformizer := f;
                Append(~evaluations, coeff_field!1);
            end if;
        end if;
    end for;
    return evaluations;
end function;

function LeadingCoefficientVectors(fs, p)
/*
Given a sequence of functions fs in a function field, and a place p on the same function field,
returns nth coefficient of the Laurent series expansion of the fs at p, where -n is the
minimum valuation of all fs at p.
*/
    coefficients := LeadingCoefficients(fs, p);
    F := ConstantField(FunctionField(p));
    return Matrix([Eltseq(x, F) : x in coefficients]);
end function;

function PrecomputeLeadingCoefficients(fs, ps)
/*
Given a sequence of functions fs in a function field, and a sequence of places ps on the same function field,
returns the leading coefficients of the Laurent series expansion of the fs at all the places of ps.
*/
    return [*LeadingCoefficientVectors(fs, p) : p in ps*];
end function;

function VanishingMatrix(coefficients, indices)
/*
Given the precomputed leading coefficients of the functions at all places, and a sequence of indices,
returns the matrix whose nullspace represents the functions vanishing at all places of indices.
*/
    return HorizontalJoin(<coefficients[i] : i in indices>);
end function;

function RiemannRochBasis(fs, coefficients, indices)
/*
Given a sequence of functions fs in a function field, the sequence coefficients returned by
PrecomputeLeadingCoefficients, and a sequence of indices, returns a basis of the space of functions of fs
vanishing at all places in indices.

TODO: is there a way to speed up the last line?
*/
    matrix := VanishingMatrix(coefficients, indices);
    nullspace := NullspaceMatrix(matrix);
    return [&+[fs[j]*nullspace[i][j] : j in [1..Ncols(nullspace)]] : i in [1..Nrows(nullspace)]];
end function;

// From here, all functions take (at most) three common arguments:
//   F          a function field
//   cg_data    F`classgroup_data
//   log        F`classgroup_data`log

procedure CachePlaceData(~cg_data, place, index)
    if IsFinite(place) then
        min := Minimum(place);
        if IsDefined(cg_data`factor_basis_minima, min) then
            Append(~cg_data`factor_basis_minima[min], index);
        else
            cg_data`factor_basis_minima[min] := [index];
        end if;
    else
        Append(~cg_data`infinite_places, index);
    end if;
end procedure;

procedure AddSupportOfDivisorToFactorBasis(~cg_data, D)
    for p in Support(D) do
        if not p in cg_data`factor_basis then
            Include(~cg_data`factor_basis, p);
            for _ in [#cg_data`factor_basis_count+1..Degree(p)] do
                Append(~cg_data`factor_basis_count, 0);
            end for;
            cg_data`factor_basis_count[Degree(p)] +:= 1;

            CachePlaceData(~cg_data, p, #cg_data`factor_basis);
        end if;
    end for;
end procedure;

procedure CreateFactorBasis(F, ~cg_data, degree_bound)
/*
Creates the factor basis consisting of all places of degree at most degree_bound and a divisor of degree 1

TODO: need to add the support of the base divisor to the factor basis
*/
    places := [Places(F, i) : i in [1..degree_bound]];
    cg_data`factor_basis := {@ p : p in x, x in places @};
    cg_data`factor_basis_count := [#x : x in places];

    cg_data`factor_basis_minima := AssociativeArray();
    cg_data`infinite_places := [];
    for i -> p in cg_data`factor_basis do
        CachePlaceData(~cg_data, p, i);
    end for;

    degree_one_divisor := DivisorOfDegreeOne(F);
    AddSupportOfDivisorToFactorBasis(~cg_data, degree_one_divisor);
end procedure;

procedure InitializeClassGroupData(F)
    F`classgroup_data := rec<classgroup_data_format |
        is_complete := false,
        log := rec<log_format | relation_time := 0, rank_time := 0, h_time := 0, functions_tried := 0, relations_found := 0>
    >;
    
    // Compute best size for factor basis by estimating the number of functions needed to determine all relations
    fb_bound := ClassGroupGenerationBound(F);
    vprintf ClassGroup: "Counting places of degree at most %o...\n", fb_bound;
    P := [NumberOfPlacesDegECF(F, i) : i in [1..fb_bound]] cat [#ConstantField(F)^d div d : d in [fb_bound+1..Genus(F)]];
    estimated_cost := [&+[P[j] : j in [1..i]] * NumberOfSmoothDivisors(Genus(F), Genus(F), P) / NumberOfSmoothDivisors(Genus(F), i, P) : i in [fb_bound..Genus(F)]];
    _, minimum_index := Minimum(estimated_cost);
    fb_bound +:= minimum_index - 1;

    vprintf ClassGroup: "Creating factor basis of places of degree at most %o...\n", fb_bound;
    CreateFactorBasis(F, ~F`classgroup_data, fb_bound);

    F`classgroup_data`relations := SparseMatrix(0, #F`classgroup_data`factor_basis);
    F`classgroup_data`empty_columns := SequenceToSet([1..#F`classgroup_data`factor_basis]);
    F`classgroup_data`rank := 0;

    h_approx := ClassNumberApproximation(F, Sqrt(2) - 1.001);
    h_lower_bound := Ceiling(h_approx / (Sqrt(2) - 0.001));
    h_upper_bound := Floor(h_approx * (Sqrt(2) - 0.001));
    F`classgroup_data`h_bounds := [h_lower_bound, h_upper_bound];
end procedure;

function GenerateRandomDivisor(cg_data, max_degree)
/*
Returns a random maximal set of places of the factor basis of total degree less than max_degree. If the relation matrix
has non-empty columns, guarantees that one of the places returned corresponds to such a column.
*/
    indices := [];
    degree := 0;

    // Add one of places corresponding to an empty column of the relation matrix
    if #cg_data`empty_columns gt 0 then
        j := Random(cg_data`empty_columns);
        Append(~indices, j);
        degree +:= Degree(cg_data`factor_basis[j]);
    end if;

    // Complete the remainder of the divisor with random places of F until the degree is as large as allowed
    while degree lt max_degree do
        max_possible_index := &+[n : i -> n in cg_data`factor_basis_count | i le max_degree - degree];
        if max_possible_index gt #indices then
            repeat
                j := Random(1, max_possible_index);
            until not j in indices;
        else
            possible_indices := [n : n in [1..max_possible_index] | not n in indices];
            if #possible_indices eq 0 then
                break;
            end if;
            j := Random(possible_indices);
        end if;
        Append(~indices, j);
        degree +:= Degree(cg_data`factor_basis[j]);
    end while;

    return indices;
end function;

function IsSupportedOnFactorBasis(cg_data, f)
/*
Returns whether the divisor of the function f is supported on the factor basis, and if so, the decomposition of the
divisor of f.
*/
    min, denom := Minimum(f, MaximalOrderFinite(Parent(f)));
    min_decomp := {x[1] : x in Factorization(min)};
    denom_decomp := {x[1] : x in Factorization(denom)};
    if not min_decomp subset Keys(cg_data`factor_basis_minima) then   // No need to test denom_decomp as we assume the pole divisor of f is supported on ps
        return false, _;
    end if;

    degree := 0;
    decomp := SparseMatrix(1, #cg_data`factor_basis);
    support := cg_data`infinite_places;
    for p in min_decomp join denom_decomp do
        support cat:= cg_data`factor_basis_minima[p];
    end for;
    
    for i in support do
        p := cg_data`factor_basis[i];
        valuation := Valuation(f, p);
        degree +:= valuation*Degree(p);
        decomp[1, i] := valuation;
    end for;
    
    if degree eq 0 then
        return true, decomp;
    else
        return false, _;
    end if;
end function;

function FindSmoothFunction(cg_data, basis)
/*
Given a sequence of functions basis, returns whether some function in the vector space spanned by the functions of basis
is supported on the factor basis. If so, also returns the decomposition of the divisor of this function.
*/
    functions_tried := 0;
    for j in [1..#basis] do
        f := basis[j];
        mult := [ConstantField(Parent(f))!0 : _ in [1..j-1]];
        while true do
            functions_tried +:= 1;
            b, rel := IsSupportedOnFactorBasis(cg_data, f);
            if b then
                return true, rel, functions_tried;
            end if;

            finished := true;
            for k in [1..#mult] do
                f +:= basis[k];
                mult[k] +:= 1;
                if not IsZero(mult[k]) then
                    finished := false;
                    break;
                end if;
            end for;
            if finished then
                break;
            end if;
        end while;
    end for;
    return false, _, functions_tried;
end function;

procedure CheckRelationLatticeCompleteness(~cg_data, ~log)
/*
Checks whether the computed relations span the whole relation lattice for the factor basis.
*/
    if cg_data`rank ne #cg_data`factor_basis - 1 then
        t := Cputime();
        cg_data`rank := Rank(cg_data`relations);
        log`latest_rank_time := Cputime(t);
        log`rank_time +:= log`latest_rank_time;
    end if;

    if cg_data`rank eq #cg_data`factor_basis - 1 then
        t := Cputime();
        cg_data`elementary_divisors := [n : n in ElementaryDivisors(cg_data`relations) | n ne 1];
        cg_data`h := &*cg_data`elementary_divisors;
        log`latest_h_time := Cputime(t);
        log`h_time +:= log`latest_h_time;

        if cg_data`h lt cg_data`h_bounds[1] then
            error "Approximation does not match up with computed order!";
        end if;
        if cg_data`h le cg_data`h_bounds[2] then
            cg_data`is_complete := true;
        end if;
    end if;
end procedure;

procedure FindRelations(~cg_data, ~log, fs)
/*
Computes relations between places of the factor basis by looking at random functions in the space spanned by the
sequence of functions fs.
*/
    vprint ClassGroup: "Precomputing function expansions...";
    coefficients := PrecomputeLeadingCoefficients(fs, cg_data`factor_basis);

    vprint ClassGroup: "Finding relations...";
    next_rank_check := log`relations_found + #cg_data`factor_basis - cg_data`rank;
    log`latest_relation_time := 0;
    log`latest_functions_tried := 0;
    log`latest_relations_found := 0;
    last_relation_start := Cputime();
    last_progress_report := Realtime();
    while true do
        indices := GenerateRandomDivisor(cg_data, #fs - 1);
        basis := RiemannRochBasis(fs, coefficients, indices);
        b, rel, functions_tried_for_basis := FindSmoothFunction(cg_data, basis);
        log`latest_functions_tried +:= functions_tried_for_basis;
        log`functions_tried +:= functions_tried_for_basis;

        if b then
            cg_data`relations := VerticalJoin(cg_data`relations, rel);
            for t in Eltseq(rel) do
                Exclude(~cg_data`empty_columns, t[2]);
            end for;
            log`latest_relations_found +:= 1;
            log`relations_found +:= 1;
        end if;

        // Logging
        if Realtime(last_progress_report) gt 10 then
            vprintf ClassGroup: "Relations found: %o/%o, smoothness probability: %o\n", log`relations_found, next_rank_check, RealField(5)!log`latest_relations_found / log`latest_functions_tried;
            last_progress_report := Realtime();
        end if;

        if log`relations_found ge next_rank_check then
            log`latest_relation_time +:= Cputime(last_relation_start);
            log`relation_time +:= Cputime(last_relation_start);

            CheckRelationLatticeCompleteness(~cg_data, ~log);

            if cg_data`is_complete then
                return;
            end if;

            relation_rate := log`latest_relations_found / log`latest_relation_time;
            if cg_data`rank lt #cg_data`factor_basis - 1 then
                next_rank_check +:= Max(
                    #cg_data`factor_basis - 1 - cg_data`rank,
                    Ceiling(5 * relation_rate * log`latest_rank_time)); // Spend at least 5 times as long searching for relations than computing rank
                vprintf ClassGroup: "Rank: %o/%o. Next check in %o relations.\n", cg_data`rank, #cg_data`factor_basis - 1, next_rank_check - log`relations_found;
            else
                next_rank_check +:= Max(1, Ceiling(5 * relation_rate * log`latest_h_time)); // Spend at least 5 times as long searching for relations than computing order
                vprintf ClassGroup: "Full rank, order: %o. Next check in %o relations.\n", cg_data`h, next_rank_check - log`relations_found;
            end if;

            last_relation_start := Cputime();
        end if;
    end while;
end procedure;

procedure ComputeClassGroupData(F : BaseDivisor := false)
/*
Computes a factor basis and relation matrix for the function field F. The function field F must be given as an extension
of the rational function field k(t) and defined over its exact constant field.
*/
    if not assigned F`classgroup_data then
        InitializeClassGroupData(F);
    end if;

    if F`classgroup_data`is_complete then
        return;
    end if;

    x := F ! BaseRing(F).1;
    if BaseDivisor cmpeq false then
        BaseDivisor := DivisorGroup(F) ! 0;
        m := Ceiling((2 * Genus(F) + #F`classgroup_data`factor_basis_count) / Degree(x));
    else
        m := 0;
    end if;

    AddSupportOfDivisorToFactorBasis(~F`classgroup_data, BaseDivisor + m * Denominator(Divisor(x)));
    if #F`classgroup_data`factor_basis gt Ncols(F`classgroup_data`relations) then
        F`classgroup_data`empty_columns join:= SequenceToSet([Ncols(F`classgroup_data`relations)+1 .. #F`classgroup_data`factor_basis]);
        F`classgroup_data`relations := HorizontalJoin(
            F`classgroup_data`relations,
            SparseMatrix(Nrows(F`classgroup_data`relations), #F`classgroup_data`factor_basis - Ncols(F`classgroup_data`relations)));
    end if;

    basis, basis_mult := ShortBasis(BaseDivisor);
    fs := [basis[i] * x^j : j in [0..basis_mult[i] + m], i in [1..#basis]];

    expected_degree := Degree(BaseDivisor) + m * Degree(x) - #fs + 1;
    P := F`classgroup_data`factor_basis_count cat [#ConstantField(F)^d div d : d in [#F`classgroup_data`factor_basis_count+1..expected_degree]];
    expected_smoothness_probability := NumberOfSmoothDivisors(expected_degree, #F`classgroup_data`factor_basis_count, P) / NumberOfSmoothDivisors(expected_degree, expected_degree, P);
    vprintf ClassGroup: "Base divisor has degree %o and Riemann-Roch dimension %o\nLooking for %o-smooth divisors of degree %o\nExpected smoothness probability: %o\n",
        (Degree(BaseDivisor) + m * Degree(x)), #fs, #F`classgroup_data`factor_basis_count, expected_degree, RealField(5)!expected_smoothness_probability;
    
    FindRelations(~F`classgroup_data, ~F`classgroup_data`log, fs);
end procedure;

intrinsic ClassGroup(F::FldFun : BaseDivisor := false) -> GrpAb
{The divisor class group of the function field F}
    F := ConstantFieldExtension(F, ExactConstantField(F)); // Ensure that the constant field of F is equal to its exact constant field
    F := RationalExtensionRepresentation(F);
    ComputeClassGroupData(F : BaseDivisor := BaseDivisor);
    
    return AbelianGroup(F`classgroup_data`elementary_divisors);
end intrinsic;
