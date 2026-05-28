procedure TimeHasFunctionOfDegreeAtMost(FF, d, label, output_file : HessMaximumTime := 100, StopAfterFirst := true)
    q := #ExactConstantField(FF);
    g := Genus(FF);

    _, timing_data := CAHasFunctionOfDegreeAtMost(FF, d : TimingData := true, StopAfterFirst := StopAfterFirst);
    fprintf output_file, "\"%o\",%o,%o,%o,Linear algebra,%o,%o,%o,%o,%o,%o,%o\n", label, q, g, d,
        timing_data`place_degree_bound, timing_data`places, timing_data`divisors, timing_data`place_enumeration_time,
        timing_data`expansions_time, timing_data`riemann_roch_time, timing_data`timeout;

    _, timing_data := CAHasFunctionOfDegreeAtMost(FF, d : Method := "Hess", MaximumTime := HessMaximumTime, TimingData := true, StopAfterFirst := StopAfterFirst);
    fprintf output_file, "\"%o\",%o,%o,%o,Hess,%o,%o,%o,%o,%o,%o,%o\n", label, q, g, d,
        timing_data`place_degree_bound, timing_data`places, timing_data`divisors, timing_data`place_enumeration_time,
        timing_data`expansions_time, timing_data`riemann_roch_time, timing_data`timeout;
end procedure;

procedure TimeClassGroup(FF, label, output_file : MaximumTime := Infinity())
    q := #ExactConstantField(FF);
    g := Genus(FF);

    start_time := Cputime();
    G := CAClassGroup(FF : MaximumTime := MaximumTime);
    if G cmpeq -1 then
        // Timeout
        fprintf output_file, "\"%o\",%o,%o,Linear algebra,,%o,true\n", label, q, g, Cputime(start_time);
    else
        fprintf output_file, "\"%o\",%o,%o,Linear algebra,%o,%o,false\n", label, q, g, Order(TorsionSubgroup(G)), Cputime(start_time);
    end if;

    // As the ClassNumber intrinsic does not provide a way to abort the calculation, we instead run the computation in
    // a separate process, which is killed after the provided maximum time.
    
    // Set up a socket for inter-process communication
    server_socket := Socket( : LocalHost := "localhost");
    t := SocketInformation(server_socket);
    host := t[1];
    port := t[2];

    pid := Fork();

    if pid eq 0 then
        // Child process
        client_socket := Socket(host, port);

        if IsFinite(MaximumTime) then
            Alarm(MaximumTime);
        end if;

        h := ClassNumber(FF);
        Write(client_socket, IntegerToString(h));
        quit;
    else
        // Parent process
        C := WaitForConnection(server_socket);
        start_time := Realtime();
        b, msg := ReadCheck(C);
        if not b or IsEof(msg) then
            // Child process terminated before finishing
            fprintf output_file, "\"%o\",%o,%o,Magma,,%o,true\n", label, q, g, Realtime(start_time);
        else
            fprintf output_file, "\"%o\",%o,%o,Magma,%o,%o,false\n", label, q, g, msg, Realtime(start_time);
        end if;
    end if;
end procedure;
