%generate wiring information for 32ch tetrode probe from Rice 
%standard format for columns in probewiring:
%column 1: channel number
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_Rice32Tinv

probewiring=[
    4	0	0	0	1
6	-14.2800000000000	0	14.2800000000000	1
17	14.2800000000000	0	14.2800000000000	1
19	0	0	28.5600000000000	1
1	0	0	150	1
2	-14.2800000000000	0	164.280000000000	1
21	14.2800000000000	0	164.280000000000	1
23	0	0	178.560000000000	1
10	0	0	300	1
8	-14.2800000000000	0	314.280000000000	1
30	14.2800000000000	0	314.280000000000	1
29	0	0	328.560000000000	1
14	0	0	450	1
12	-14.2800000000000	0	464.280000000000	1
27	14.2800000000000	0	464.280000000000	1
25	0	0	478.560000000000	1
5	0	0	600	1
7	-14.2800000000000	0	614.280000000000	1
16	14.2800000000000	0	614.280000000000	1
18	0	0	628.560000000000	1
0	0	0	750	1
3	-14.2800000000000	0	764.280000000000	1
20	14.2800000000000	0	764.280000000000	1
22	0	0	778.560000000000	1
11	0	0	900	1
9	-14.2800000000000	0	914.280000000000	1
31	14.2800000000000	0	914.280000000000	1
28	0	0	928.560000000000	1
15	0	0	1050	1
13	-14.2800000000000	0	1064.28000000000	1
26	14.2800000000000	0	1064.28000000000	1
24	0	0	1078.56000000000	1];