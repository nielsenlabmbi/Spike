%first column = channel number on probe.  Convention is ch. 1 is at top left of PCB
%and number wrap around in counter-clockwise fashion.

%second column = x coordinate, horizontal distance to reference electrode in microns.  Convention is
%for reference electrode to be at bottom left corner of array when array
%points up.

%third column = y coordinate.  Convention is for this coordinate to be zero
%unless the probe is 3D.

%forth column = z coordinate, vertical distance to reference electrode in microns.  Convention is
%for reference electrode to be at bottom left corner of array when array
%points up.

%fifth column = shaft number.  Convention is for shaft 1 to be the
%right-most shaft when array points up.

probewiring=[
1	0	0	600	4
2	0	0	550	4
3	0	0	500	4
4	0	0	450	4
5	0	0	400	4
6	0	0	350	4
7	0	0	300	4
8	0	0	250	4
9	0	0	200	4
10	0	0	150	4
11	0	0	100	4
12	0	0	50	4
13	0	0	0	4
14	0	0	650	4
15	0	0	700	4
16	3.9	0	750	4
17	10	0	775	4
18	20	0	725	4
19	20	0	675	4
20	20	0	25	4
21	20	0	75	4
22	20	0	125	4
23	20	0	175	4
24	20	0	225	4
25	20	0	275	4
26	20	0	325	4
27	20	0	375	4
28	20	0	425	4
29	20	0	475	4
30	20	0	525	4
31	20	0	575	4
32	20	0	625	4
33	150	0	100	3
34	150	0	50	3
35	150	0	0	3
36	150	0	-50	3
37	150	0	-100	3
38	150	0	-150	3
39	150	0	-200	3
40	150	0	-250	3
41	150	0	-300	3
42	150	0	-350	3
43	150	0	-400	3
44	150	0	-450	3
45	150	0	-500	3
46	150	0	150	3
47	150	0	200	3
48	153.9	0	250	3
49	160	0	275	3
50	170	0	225	3
51	170	0	175	3
52	170	0	-475	3
53	170	0	-425	3
54	170	0	-375	3
55	170	0	-325	3
56	170	0	-275	3
57	170	0	-225	3
58	170	0	-175	3
59	170	0	-125	3
60	170	0	-75	3
61	170	0	-25	3
62	170	0	25	3
63	170	0	75	3
64	170	0	125	3
65	300	0	-400	2
66	300	0	-450	2
67	300	0	-500	2
68	300	0	-550	2
69	300	0	-600	2
70	300	0	-650	2
71	300	0	-700	2
72	300	0	-750	2
73	300	0	-800	2
74	300	0	-850	2
75	300	0	-900	2
76	300	0	-950	2
77	300	0	-1000	2
78	300	0	-350	2
79	300	0	-300	2
80	303.9	0	-250	2
81	310	0	-225	2
82	320	0	-275	2
83	320	0	-325	2
84	320	0	-975	2
85	320	0	-925	2
86	320	0	-875	2
87	320	0	-825	2
88	320	0	-775	2
89	320	0	-725	2
90	320	0	-675	2
91	320	0	-625	2
92	320	0	-575	2
93	320	0	-525	2
94	320	0	-475	2
95	320	0	-425	2
96	320	0	-375	2
97	450	0	-900	1
98	450	0	-950	1
99	450	0	-1000	1
100	450	0	-1050	1
101	450	0	-1100	1
102	450	0	-1150	1
103	450	0	-1200	1
104	450	0	-1250	1
105	450	0	-1300	1
106	450	0	-1350	1
107	450	0	-1400	1
108	450	0	-1450	1
109	450	0	-1500	1
110	450	0	-850	1
111	450	0	-800	1
112	453.9	0	-750	1
113	460	0	-725	1
114	470	0	-775	1
115	470	0	-825	1
116	470	0	-1475	1
117	470	0	-1425	1
118	470	0	-1375	1
119	470	0	-1325	1
120	470	0	-1275	1
121	470	0	-1225	1
122	470	0	-1175	1
123	470	0	-1125	1
124	470	0	-1075	1
125	470	0	-1025	1
126	470	0	-975	1
127	470	0	-925	1
128	470	0	-875	1

];


tipelectrode=30;    %nearest tip-electrode vertical distance in microns.

connector_position='top';

find_probewiring



s=[];
s.channels=probewiring(:,1);
s.x=probewiring(:,2);   %the reference electrode is always the top right channel when the probes are pointing up.
s.y=probewiring(:,3);
s.z=probewiring(:,4); s.z=s.z-min(s.z);
s.shaft=probewiring(:,5);
s.tipelectrode=tipelectrode;

if strcmp(headstage_source,'Intan')
    s.channels=s.channels-1;
end

%To plot the labeled channels:
figure(2)
clf
plot(s.x,s.z,'sqr', 'MarkerSize',11)
hold on
for i=2:2:size(probewiring,1)
text(s.x(i)-5,s.z(i),num2str(s.channels(i)),'FontSize',10)
end
axis([min(s.x)-50 max(s.x)+50 min(s.z)-50 max(s.z)+50])
axis equal
set(gca,'FontSize',10,'TickDir','out')
