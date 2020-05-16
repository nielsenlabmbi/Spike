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
1	0	0	1150	2
2	0	0	1100	2
3	0	0	1050	2
4	0	0	1000	2
5	0	0	950	2
6	0	0	900	2
7	0	0	850	2
8	0	0	800	2
9	0	0	750	2
10	0	0	700	2
11	0	0	650	2
12	0	0	600	2
13	0	0	550	2
14	0	0	500	2
15	0	0	450	2
16	0	0	400	2
17	0	0	350	2
18	0	0	300	2
19	0	0	250	2
20	0	0	200	2
21	0	0	150	2
22	0	0	100	2
23	0	0	50	2
24	0	0	0	2
25	0	0	1200	2
26	0	0	1250	2
27	0	0	1300	2
28	0	0	1350	2
29	0	0	1400	2
30	0	0	1450	2
31	0	0	1500	2
32	0	0	1550	2
33	25	0	1575	2
34	25	0	1525	2
35	25	0	1475	2
36	25	0	1425	2
37	25	0	1375	2
38	25	0	1325	2
39	25	0	1275	2
40	25	0	1225	2
41	25	0	1175	2
42	25	0	25	2
43	25	0	75	2
44	25	0	125	2
45	25	0	175	2
46	25	0	225	2
47	25	0	275	2
48	25	0	325	2
49	25	0	375	2
50	25	0	425	2
51	25	0	475	2
52	25	0	525	2
53	25	0	575	2
54	25	0	625	2
55	25	0	675	2
56	25	0	725	2
57	25	0	775	2
58	25	0	825	2
59	25	0	875	2
60	25	0	925	2
61	25	0	975	2
62	25	0	1025	2
63	25	0	1075	2
64	25	0	1125	2
65	500	0	1150	1
66	500	0	1100	1
67	500	0	1050	1
68	500	0	1000	1
69	500	0	950	1
70	500	0	900	1
71	500	0	850	1
72	500	0	800	1
73	500	0	750	1
74	500	0	700	1
75	500	0	650	1
76	500	0	600	1
77	500	0	550	1
78	500	0	500	1
79	500	0	450	1
80	500	0	400	1
81	500	0	350	1
82	500	0	300	1
83	500	0	250	1
84	500	0	200	1
85	500	0	150	1
86	500	0	100	1
87	500	0	50	1
88	500	0	0	1
89	500	0	1200	1
90	500	0	1250	1
91	500	0	1300	1
92	500	0	1350	1
93	500	0	1400	1
94	500	0	1450	1
95	500	0	1500	1
96	500	0	1550	1
97	525	0	1575	1
98	525	0	1525	1
99	525	0	1475	1
100	525	0	1425	1
101	525	0	1375	1
102	525	0	1325	1
103	525	0	1275	1
104	525	0	1225	1
105	525	0	1175	1
106	525	0	25	1
107	525	0	75	1
108	525	0	125	1
109	525	0	175	1
110	525	0	225	1
111	525	0	275	1
112	525	0	325	1
113	525	0	375	1
114	525	0	425	1
115	525	0	475	1
116	525	0	525	1
117	525	0	575	1
118	525	0	625	1
119	525	0	675	1
120	525	0	725	1
121	525	0	775	1
122	525	0	825	1
123	525	0	875	1
124	525	0	925	1
125	525	0	975	1
126	525	0	1025	1
127	525	0	1075	1
128	525	0	1125	1

];


tipelectrode=65.1;    %nearest tip-electrode vertical distance in microns.

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
for i=1:size(probewiring,1)
text(s.x(i)-5,s.z(i),num2str(s.channels(i)),'FontSize',9)
end
axis([min(s.x)-50 max(s.x)+50 min(s.z)-50 max(s.z)+50])
axis equal
set(gca,'FontSize',10,'TickDir','out')
