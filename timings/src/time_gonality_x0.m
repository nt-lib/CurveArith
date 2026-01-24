output_file := Open("output/gonality_x0.csv", "w");

// Format is: <level, q, d>
queue := [
    <38, 5, 3>, <44, 5, 3>, <53, 7, 3>, <61, 3, 3>, <76, 5, 5>, <82, 5, 5>, <84, 5, 5>, <86, 3, 5>, <93, 5, 5>,
    <99, 5, 5>, <102, 5, 7>, <106, 7, 7>, <108, 5, 5>, <109, 3, 4>, <112, 3, 5>, <113, 3, 5>, <114, 5, 7>, <115, 3, 5>,
    <116, 3, 5>, <117, 5, 5>, <118, 3, 5>, <122, 3, 5>, <127, 3, 5>, <128, 3, 5>, <130, 3, 7>, <132, 5, 7>, <134, 3, 7>,
    <136, 5, 7>, <137, 3, 5>, <140, 3, 7>, <144, 5, 5>, <147, 5, 5>, <148, 5, 7>, <150, 7, 7>, <151, 5, 5>, <152, 3, 7>,
    <153, 5, 7>, <154, 3, 7>, <157, 3, 7>, <160, 7, 7>, <162, 5, 5>, <163, 3, 6>, <169, 5, 5>, <170, 3, 7>, <172, 3, 7>,
    <175, 2, 7>, <176, 3, 7>, <178, 3, 7>, <179, 5, 5>, <180, 7, 6>, <181, 3, 5>, <187, 2, 7>, <189, 2, 7>, <192, 5, 7>,
    <193, 3, 7>, <196, 5, 7>, <197, 3, 7>, <198, 5, 7>, <200, 3, 7>, <201, 2, 7>, <217, 2, 7>, <229, 3, 7>, <233, 2, 7>,
    <241, 2, 7>, <247, 2, 7>]; // Takes too much memory: <277, 5, 7>

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
            SetMemoryLimit(20 * 10^9);
            client_socket := Socket(host, port);
            n, q, d := Explode(queue[i]);
            C := ModularCurveQuotient(n, []);
            C := ChangeRing(C, GF(q));
            FF := AlgorithmicFunctionField(FunctionField(C));
            TimeHasFunctionOfDegreeAtMost(FF, d, n, output_file);
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
