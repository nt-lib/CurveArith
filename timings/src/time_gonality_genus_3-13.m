output_file := Open("output/gonality_genus_3-13.csv", "w");

// Format is: <genus, q, Ceiling(#X(F_q) / (q + 1))>
testing_parameters := [
    <4, 3, 1>, <4, 7, 1>, <4, 13, 1>, <4, 23, 1>, <4, 31, 1>, <4, 59, 1>, <4, 89, 1>,
    <6, 3, 1>, <6, 7, 1>, <6, 13, 1>, <6, 23, 1>, <6, 31, 1>,
    <8, 3, 1>, <8, 7, 1>, <8, 13, 1>,
    <10, 3, 1>, <10, 7, 1>,
    <4, 3, 2>, <4, 7, 2>, <4, 13, 2>, <4, 23, 2>, <4, 31, 2>, <4, 59, 2>, <4, 89, 2>,
    <6, 3, 2>, <6, 7, 2>, <6, 13, 2>, <6, 23, 2>, <6, 31, 2>, <6, 59, 2>, <6, 89, 2>,
    <8, 3, 2>, <8, 7, 2>, <8, 13, 2>, <8, 23, 2>, <8, 31, 2>,
    <10, 3, 2>, <10, 7, 2>,
    <11, 2, 1>, <11, 3, 1>, <11, 5, 1>,
    <12, 2, 1>, <12, 3, 1>, <12, 5, 1>,
    <13, 2, 1>, <13, 3, 1>,
    <11, 2, 2>, <11, 3, 2>, <11, 5, 2>, <11, 7, 2>,
    <12, 2, 2>, <12, 3, 2>, <12, 5, 2>, <12, 7, 2>,
    <13, 2, 2>, <13, 3, 2>, <13, 5, 2>];

queue := [x : _ in [1..5], x in testing_parameters];

// Basic multi-process implementation
server_socket := Socket( : LocalHost := "localhost");
t := SocketInformation(server_socket);
host := t[1];
port := t[2];

processes := 5;

seed := GetSeed();

i := 0;
finished := 0;
read_sockets := {};
last_progress := Realtime();
while finished lt #queue do
    for _ in [1..Minimum(processes - #read_sockets, #queue - i)] do
        i +:= 1;
        seed +:= 1;
        SetSeed(seed);
        pid := Fork();

        if pid eq 0 then
            SetMemoryLimit(5 * 10^9);
            client_socket := Socket(host, port);
            g, q, b := Explode(queue[i]);
            while true do
                try
                    C := RandomCurveByGenus(g, GF(q));
                    FF := AlgorithmicFunctionField(FunctionField(C));
                    if Ceiling(#Places(FF, 1) / (q + 1)) eq b then
                        break;
                    end if;
                catch e
                    print "Error in random curve generation";
                end try;
            end while;
            TimeHasFunctionOfDegreeAtMost(FF, (Genus(FF) + 3) div 2, &cat Split(Sprint(FF, "Magma"), "\n"), output_file : StopAfterFirst := false);
            Write(client_socket, "done");
            quit;
        else
            Include(~read_sockets, WaitForConnection(server_socket));
        end if;
    end for;

    ready := WaitForIO(Setseq(read_sockets));
    for I in ready do
        b, msg := ReadCheck(I);
        if not b or IsEof(msg) then
            "Closed without finishing!";
        end if;
        finished +:= 1;
        Exclude(~read_sockets, I);
    end for;

    if Realtime(last_progress) ge 10 then
        printf "Progress: %o/%o\n", finished, #queue;
        last_progress := Realtime();
    end if;
end while;
quit;
