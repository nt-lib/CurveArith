AttachSpec("../CurveArith.spec");

start_time := Realtime();
load "src/timing.m";
load "src/time_gonality_genus_3-13.m";
load "src/time_gonality_x0.m";
load "src/time_classgroup_genus_0-13.m";

printf "All timings finished in %o!\n", Realtime(start_time);
exit;
