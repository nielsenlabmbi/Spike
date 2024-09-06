%configuration for NeuroNexus A16x1-3mm-100-177
%(16 shanks, 100 um shank spacing 
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_A16x1

probewiring=[
    5	0	0	0	1
2	100	0	0	2
0	200	0	0	3
1	300	0	0	4
4	400	0	0	5
3	500	0	0	6
6	600	0	0	7
7	700	0	0	8
8	800	0	0	9
9	900	0	0	10
12	1000	0	0	11
11	1100	0	0	12
14	1200	0	0	13
15	1300	0	0	14
13	1400	0	0	15
10	1500	0	0	16];