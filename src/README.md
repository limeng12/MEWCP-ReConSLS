use this command to run:
./ReConSLS -inst instancefile -seed seed -cutoff cutofftime -min_bms minbms -max_bms maxbms -Waim waim

-inst:  the instance file   (necessary)
-seed:  seed of the program (necessary)
-cutoff:    the time allowed to run the program (necessary)
-min_bms:   the default setting is 8    (not necessary)
-max_bms:   the default setting is 128  (not necessary)
-Waim:      the program stops when a solution greater than or equal to the value is found (not necessary)

two types of instance file:
1. use #define format_dimacs
edge 154908 327162
e 1 2
e 1 3
e 1 4
e 1 5
...

2. use #define format_konect
edge 154908 327162
e 1 2 4
e 1 3 5
e 1 4 6
e 1 5 7
...
