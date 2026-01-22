procedure TimeHasFunctionOfDegreeAtMost(FF, d, label, output_file : HessMaximumTime := 100, StopAfterFirst := true)
    q := #ExactConstantField(FF);
    g := Genus(FF);

    _, timing_data := HasFunctionOfDegreeAtMost(FF, d : TimingData := true, StopAfterFirst := StopAfterFirst);
    fprintf output_file, "%o;%o;%o;%o;Linear algebra;%o;%o;%o;%o;%o;%o;%o\n", label, q, g, d,
        timing_data`place_degree_bound, timing_data`places, timing_data`divisors, timing_data`place_enumeration_time,
        timing_data`expansions_time, timing_data`riemann_roch_time, timing_data`timeout;

    _, timing_data := HasFunctionOfDegreeAtMost(FF, d : Method := "Hess", MaximumTime := HessMaximumTime, TimingData := true, StopAfterFirst := StopAfterFirst);
    fprintf output_file, "%o;%o;%o;%o;Hess;%o;%o;%o;%o;%o;%o;%o\n", label, q, g, d,
        timing_data`place_degree_bound, timing_data`places, timing_data`divisors, timing_data`place_enumeration_time,
        timing_data`expansions_time, timing_data`riemann_roch_time, timing_data`timeout;
end procedure;

procedure TimeClassGroup(FF, label, output_file)
    q := #ExactConstantField(FF);
    g := Genus(FF);

    fprintf output_file, "%o;%o;%o\n", label, q, g; // Allows us to know if a timeout happened, since there's no way to implement a timeout for the original class group code without killing the process

    start_time := Cputime();
    G := ClassGroup(FF);
    fprintf output_file, "%o;%o;%o;Linear algebra;%o;%o\n", label, q, g, Order(TorsionSubgroup(G)), Cputime(start_time);

    start_time := Cputime();
    h := ClassNumber(FF);
    fprintf output_file, "%o;%o;%o;Original;%o;%o\n", label, q, g, h, Cputime(start_time);
end procedure;
