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
