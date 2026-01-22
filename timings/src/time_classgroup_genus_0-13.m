output_file := Open("output/classgroup_genus_0-13.csv", "a");

// Format is: <genus, q>
testing_parameters := [<g, q> : q in [2, 3, 5, 7, 11, 13], g in [0..13]];

queue := [x : _ in [1..5], x in testing_parameters];

// Basic multi-process implementation
server_socket := Socket( : LocalHost := "localhost");
t := SocketInformation(server_socket);
host := t[1];
port := t[2];

processes := 5;
timeout := 1000;

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
            g, q := Explode(queue[i]);
            while true do
                try
                    C := RandomCurveByGenus(g, GF(q));
                    FF := AlgorithmicFunctionField(FunctionField(C));
                    break;
                catch e
                    printf "Error in random curve generation (g=%o, q=%o)\n", g, q;
                end try;
            end while;
            Alarm(timeout);
            TimeClassGroup(FF, &cat Split(Sprint(FF, "Magma"), "\n"), output_file);
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
