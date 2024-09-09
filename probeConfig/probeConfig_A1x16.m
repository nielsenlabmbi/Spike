%configuration for NeuroNexus A1x16-3mm-100-177
%1 shanks, 100 um shank spacing 
%note - this wiring does not take the mapping between connector and
%headstage into account
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_A1x16

probewiring=[
    5	0	0	0	1
10	0	0	-100	1
2	0	0	-200	1
13	0	0	-300	1
0	0	0	-400	1
15	0	0	-500	1
1	0	0	-600	1
14	0	0	-700	1
4	0	0	-800	1
11	0	0	-900	1
3	0	0	-1000	1
12	0	0	-1100	1
6	0	0	-1200	1
9	0	0	-1300	1
7	0	0	-1400	1
8	0	0	-1500	1];