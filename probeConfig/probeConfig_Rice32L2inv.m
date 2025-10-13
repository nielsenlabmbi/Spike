%generate wiring information for 32ch linear probe from Rice 
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_Rice32L2inv

probewiring=[
24	25	0	0	1
11	0	0	25	1
26	25	0	50	1
9	0	0	75	1
17	25	0	100	1
6	0	0	125	1
19	25	0	150	1
4	0	0	175	1
21	25	0	200	1
2	0	0	225	1
23	25	0	250	1
1	0	0	275	1
31	25	0	300	1
8	0	0	325	1
29	25	0	350	1
10	0	0	375	1
27	25	0	400	1
12	0	0	425	1
25	25	0	450	1
14	0	0	475	1
16	25	0	500	1
7	0	0	525	1
18	25	0	550	1
5	0	0	575	1
20	25	0	600	1
3	0	0	625	1
22	25	0	650	1
0	0	0	675	1
30	25	0	700	1
15	0	0	725	1
28	25	0	750 1
13	0	0	775	1];