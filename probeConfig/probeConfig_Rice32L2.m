%generate wiring information for 32ch linear probe from Rice 
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_Rice32L2

probewiring=[
    8	25	0	0	1
27	0	0	25	1
10	25	0	50	1
25	0	0	75	1
1	25	0	100	1
22	0	0	125	1
3	25	0	150	1
20	0	0	175	1
5	25	0	200	1
18	0	0	225	1
7	25	0	250	1
17	0	0	275	1
15	25	0	300	1
24	0	0	325	1
13	25	0	350	1
26	0	0	375	1
11	25	0	400	1
28	0	0	425	1
9	25	0	450	1
30	0	0	475	1
0	25	0	500	1
23	0	0	525	1
2	25	0	550	1
21	0	0	575	1
4	25	0	600	1
19	0	0	625	1
6	25	0	650	1
16	0	0	675	1
14	25	0	700	1
31	0	0	725	1
12	25	0	750 1
29	0	0	775	1];