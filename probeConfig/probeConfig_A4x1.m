%configuration for NeuroNexus A4x1-tet-3mm-150-121
%4 shanks, 1 tetrode each shank, 150 um shank spacing
%assumption is 25um spacing between sites
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_A4x1

probewiring=[
    0	0	0	0	1
5	-18	0	18	1
1	18	0	18	1
2	0	0	36	1
6	150	0	0	2
5	132	0	18	2
7	168	0	18	2
3	150	0	36	2
12	300	0	0	3
8	282	0	18	3
11	318	0	18	3
9	300	0	36	3
13	450	0	0	4
14	432	0	18	4
10	468	0	18	4
15	450	0	36	4];