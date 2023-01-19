%generate wiring information for 32ch linear probe from Rice 
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_Rice32L

probewiring=[
    9	0	0	0	1
30	0	0	60	1
11	0	0	120	1
28	0	0	180	1
13	0	0	240	1
26	0	0	300	1
14	0	0	360	1
24	0	0	420	1
7	0	0	480	1
17	0	0	540	1
5	0	0	600	1
18	0	0	660	1
3	0	0	720	1
20	0	0	780	1
1	0	0	840	1
22	0	0	900	1
8	0	0	960	1
31	0	0	1020	1
10	0	0	1080	1
29	0	0	1140	1
12	0	0	1200	1
27	0	0	1260	1
15	0	0	1320	1
25	0	0	1380	1
6	0	0	1440	1
16	0	0	1500	1
4	0	0	1560	1
19	0	0	1620	1
2	0	0	1680	1
21	0	0	1740	1
0	0	0	1800	1
23	0	0	1860	1];