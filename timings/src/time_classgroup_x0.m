output_file := Open("output/classgroup_x0.csv", "w");

// Format is: <level, q>
queue := [<n, q> : q in [2, 3, 7, 17, 31, 59, 97], n in [1..150] | not n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25] and not n mod q eq 0];

// Basic multi-process implementation
server_socket := Socket( : LocalHost := "localhost");
t := SocketInformation(server_socket);
host := t[1];
port := t[2];

processes := 10;
timeout := 3600;

i := 0;
finished := 0;
read_sockets := {};
last_progress := Realtime();
while finished lt 2 * #queue do
    for _ in [1..Minimum(processes - #read_sockets div 2, #queue - i)] do
        i +:= 1;
        pid := Fork();

        if pid eq 0 then
            SetMemoryLimit(2 * 10^9);
            n, q := Explode(queue[i]);
            try
                C := ModularCurveQuotient(n, []);
                C := ChangeRing(C, GF(q));
                FF := AlgorithmicFunctionField(FunctionField(C));
            catch e
                printf "Error in curve generation: %o", e;
                Socket(host, port);
                Socket(host, port);
                quit;
            end try;

            pid := Fork();
            if pid eq 0 then
                client_socket := Socket(host, port);
                Alarm(timeout);
                TimeClassGroupLinearAlgebra(FF, n, output_file);
                Write(client_socket, "done");
                quit;
            else
                client_socket := Socket(host, port);
                Alarm(timeout);
                TimeClassGroupOriginal(FF, n, output_file);
                Write(client_socket, "done");
                quit;
            end if;
        else
            Include(~read_sockets, WaitForConnection(server_socket));
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
        printf "Progress: %o/%o\n", finished, 2 * #queue;
        last_progress := Realtime();
    end if;
end while;
quit;
