%generate wiring information for single electrode
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_PlexonSProbeVH

probewiring=[
23	0	0	750	1
22	50	0	750	1
21	0	0	700	1
20	50	0	700	1
19	0	0	650	1
18	50	0	650	1
17	0	0	600	1
16	50	0	600	1
15	0	0	550	1
14	50	0	550	1
13	0	0	500	1
12	50	0	500	1
11	0	0	450	1
10	50	0	450	1
9	0	0	400	1
8	50	0	400	1
24	0	0	350	1
25	50	0	350	1
26	0	0	300	1
27	50	0	300	1
28	0	0	250	1
29	50	0	250	1
30	0	0	200	1
31	50	0	200	1
0	0	0	150	1
1	50	0	150	1
2	0	0	100	1
3	50	0	100	1
4	0	0	50	1
5	50	0	50	1
6	0	0	0	1
7	50	0	0	1
];