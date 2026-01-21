declare verbose ClassGroup, 1;

declare attributes FldFunG: classgroup_data;

classgroup_data_format := recformat<
    factor_basis,           // Factor basis
    factor_basis_degree,    // Largest degree d such that the factor basis contains all places of degree <= d
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

function NumberOfSmoothDivisorsDistinctSupport(n, m, P)
/*
Given two integers n and m and a sequence P containing the number of places of degree <= min(n, m) on a function field F,
returns the number of effective divisors on F of degree n, consisting only of places of degree <= m and with distinct
support.
*/
    if n lt 0 or m lt 0 then
        return 0;
    end if;

    // psi[i+1][j+1] is the number of i-smooth divisors of degree j with distinct support
    psi := [[1] cat [0 : _ in [1..n]]];
    for i in [1..m] do
        row := [1];
        for j in [1..n] do
            count := 0;
            for k in [0..Minimum(P[i], j div i)] do
                count +:= Binomial(P[i], k) * psi[i][j - k*i + 1];
            end for;
            Append(~row, count);
        end for;
        Append(~psi, row);
    end for;

    return psi[m+1][n+1];
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

procedure CreateFactorBasis(F, ~cg_data, degree_bound, base_divisor)
/*
Creates the factor basis consisting of all places of degree at most degree_bound and the supports of base_divisor
and a divisor of degree one.
*/
    // Enumerate all places of degree at most degree_bound
    places := [Places(F, i) : i in [1..degree_bound]];
    cg_data`factor_basis := [p : p in x, x in places];
    cg_data`factor_basis_degree := degree_bound;
    cg_data`factor_basis_count := [#x : x in places];

    // Add supports of base_divisor and a divisor of degree one
    degree_one_divisor := DivisorOfDegreeOne(F);
    additional_places := Setseq(Seqset(Support(degree_one_divisor) cat Support(base_divisor)));
    Sort(~additional_places, func<x, y | Degree(x) - Degree(y)>); // Factor basis must be sorted by degree
    for p in additional_places do
        if Degree(p) gt degree_bound then
            Append(~cg_data`factor_basis, p);
            for _ in [#cg_data`factor_basis_count+1 .. Degree(p)] do
                Append(~cg_data`factor_basis_count, 0);
            end for;
            cg_data`factor_basis_count[Degree(p)] +:= 1;
        end if;
    end for;

    // Cache underlying places in k[t]
    cg_data`factor_basis_minima := AssociativeArray();
    cg_data`infinite_places := [];
    for i -> p in cg_data`factor_basis do
        if IsFinite(p) then
            min := Minimum(p);
            if IsDefined(cg_data`factor_basis_minima, min) then
                Append(~cg_data`factor_basis_minima[min], i);
            else
                cg_data`factor_basis_minima[min] := [i];
            end if;
        else
            Append(~cg_data`infinite_places, i);
        end if;
    end for;

    // Reset relation matrix
    cg_data`relations := SparseMatrix(0, #cg_data`factor_basis);
    cg_data`empty_columns := SequenceToSet([1..#cg_data`factor_basis]);
    cg_data`rank := 0;
end procedure;

procedure InitializeClassGroupData(F)
/*
Initialize the classgroup_data record of the function field F.
*/
    F`classgroup_data := rec<classgroup_data_format |
        factor_basis_degree := 0,
        factor_basis_count := [],
        is_complete := false,
        log := rec<log_format | relation_time := 0, rank_time := 0, h_time := 0, functions_tried := 0, relations_found := 0>
    >;

    h_approx := ClassNumberApproximation(F, Sqrt(2) - 1.001);
    h_lower_bound := Ceiling(h_approx / (Sqrt(2) - 0.001));
    h_upper_bound := Floor(h_approx * (Sqrt(2) - 0.001));
    F`classgroup_data`h_bounds := [h_lower_bound, h_upper_bound];
end procedure;

function GenerateRandomDivisor(cg_data, min_degree, max_degree)
/*
Returns a random maximal set of places of the factor basis of total degree between min_degree and max_degree. If the
relation matrix has non-empty columns, guarantees that one of the places returned corresponds to such a column.
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
    while degree lt min_degree do
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
        // Start at a random function in the vector space, to avoid repetition if we land in the same space twice
        for k in [j+1..#basis] do
            f +:= Random(ConstantField(Parent(f))) * basis[k];
        end for;
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
        cg_data`h := #cg_data`elementary_divisors eq 0 select 1 else &*cg_data`elementary_divisors;
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

procedure FindRelations(~cg_data, ~log, fs, subtraction_degree)
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
        indices := GenerateRandomDivisor(cg_data, subtraction_degree, #fs - 1);
        basis := RiemannRochBasis(fs, coefficients, indices);
        b, rel, functions_tried_for_basis := FindSmoothFunction(cg_data, basis);
        log`latest_functions_tried +:= functions_tried_for_basis;
        log`functions_tried +:= functions_tried_for_basis;

        if b then
            // Relation found
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

            relation_rate := log`latest_relations_found / Maximum(0.01, log`latest_relation_time);
            if cg_data`rank lt #cg_data`factor_basis - 1 then
                next_rank_check +:= Max(
                    #cg_data`factor_basis - 1 - cg_data`rank,
                    Ceiling(5 * relation_rate * Maximum(0.01, log`latest_rank_time))); // Spend at least 5 times as long searching for relations than computing rank
                vprintf ClassGroup: "Rank: %o/%o. Next check in %o relations.\n", cg_data`rank, #cg_data`factor_basis - 1, next_rank_check - log`relations_found;
            else
                next_rank_check +:= Max(1, Ceiling(5 * relation_rate * Maximum(0.01, log`latest_h_time))); // Spend at least 5 times as long searching for relations than computing order
                vprintf ClassGroup: "Full rank, order: %o. Next check in %o relations.\n", cg_data`h, next_rank_check - log`relations_found;
            end if;

            last_relation_start := Cputime();
        end if;
    end while;
end procedure;

procedure ComputeClassGroupData(F : BaseDivisor := false, FactorBasisDegree := -1)
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
    
    if FactorBasisDegree lt 1 then
        FactorBasisDegree := ClassGroupGenerationBound(F);

        // Ensure that the factor basis has at least 20 places, as having a very small factor basis can lead to issues
        place_counts := [NumberOfPlacesDegECF(F, i) : i in [1..FactorBasisDegree]];
        while &+place_counts lt 20 do
            FactorBasisDegree +:= 1;
            Append(~place_counts, NumberOfPlacesDegECF(F, FactorBasisDegree));
        end while;
    end if;

    if BaseDivisor cmpeq false then
        base_divisor := Denominator(Divisor(F ! BaseRing(F).1));
    else
        base_divisor := BaseDivisor;
    end if;

    // Check if the base divisor is supported on current factor basis
    base_divisor_in_factor_basis := true;
    for p in Support(base_divisor) do
        if Degree(p) gt F`classgroup_data`factor_basis_degree then
            if Degree(p) gt #F`classgroup_data`factor_basis_count then
                base_divisor_in_factor_basis := false;
                break;
            end if;

            start_index := &+[F`classgroup_data`factor_basis_count[i] : i in [1 .. Degree(p)-1]];
            for i in [1..F`classgroup_data`factor_basis_count[Degree(p)]] do
                if F`classgroup_data`factor_basis[start_index + i] eq p then
                    continue p;
                end if;
            end for;
            base_divisor_in_factor_basis := false;
            break;
        end if;
    end for;

    // Enlarge factor basis if necessary
    if F`classgroup_data`factor_basis_degree lt FactorBasisDegree or not base_divisor_in_factor_basis then
        vprintf ClassGroup: "Creating factor basis of places of degree at most %o...\n", FactorBasisDegree;
        CreateFactorBasis(F, ~F`classgroup_data, FactorBasisDegree, base_divisor);
    end if;

    // If BaseDivisor is not user-specified, take a large enough multiple of the infinity divisor to (hopefully)
    // ensure all relations can be found
    if BaseDivisor cmpeq false then
        m := Ceiling((2*Genus(F) + #F`classgroup_data`factor_basis_count) / Degree(base_divisor));
        base_divisor := m * base_divisor;
    end if;

    fs := Basis(base_divisor);

    // Set the degree of the subtraction divisors in the relation finding procedure so that the number of relations
    // the procedure is expected to find is large enough. This is usually taken to be l(base_divisor) - 1, but can
    // be smaller when there are very few rational places.
    place_counts := [d le F`classgroup_data`factor_basis_degree select F`classgroup_data`factor_basis_count[d] else #ConstantField(F)^d div d : d in [1..Degree(base_divisor)]];
    subtraction_degree := #fs;
    expected_smooth_functions := 0;
    repeat
        subtraction_degree -:= 1;
        unknown_degree := Degree(base_divisor) - subtraction_degree;
        expected_smoothness_probability :=
            NumberOfSmoothDivisors(unknown_degree, #F`classgroup_data`factor_basis_count, F`classgroup_data`factor_basis_count)
                / NumberOfSmoothDivisors(unknown_degree, unknown_degree, place_counts);
        accessible_functions := (#ConstantField(F)^(#fs - subtraction_degree) - 1) / (#ConstantField(F) - 1)
            * NumberOfSmoothDivisorsDistinctSupport(subtraction_degree, #F`classgroup_data`factor_basis_count, F`classgroup_data`factor_basis_count);
        expected_smooth_functions +:= expected_smoothness_probability * accessible_functions;
    until expected_smooth_functions ge 10 * #F`classgroup_data`factor_basis or subtraction_degree eq #F`classgroup_data`factor_basis_count;

    vprintf ClassGroup: "Base divisor has degree %o and Riemann-Roch dimension %o\nLooking for %o-smooth divisors of degree %o\nExpected smoothness probability: %o\n",
        Degree(base_divisor), #fs, F`classgroup_data`factor_basis_degree, unknown_degree, RealField(5)!expected_smoothness_probability;
        
    FindRelations(~F`classgroup_data, ~F`classgroup_data`log, fs, subtraction_degree);
end procedure;

intrinsic ClassGroup(F::FldFun : BaseDivisor := false, FactorBasisDegree := -1) -> GrpAb
{The divisor class group of the function field F}
    F := ConstantFieldExtension(F, ExactConstantField(F)); // Ensure that the constant field of F is equal to its exact constant field
    F := RationalExtensionRepresentation(F);
    ComputeClassGroupData(F : BaseDivisor := BaseDivisor, FactorBasisDegree := FactorBasisDegree);
    
    return AbelianGroup(F`classgroup_data`elementary_divisors cat [0]);
end intrinsic;
