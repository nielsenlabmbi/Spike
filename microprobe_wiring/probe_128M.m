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
1	0	0	840	4
2	0	0	770	4
3	0	0	700	4
4	0	0	630	4
5	0	0	560	4
6	0	0	490	4
7	0	0	420	4
8	0	0	350	4
9	0	0	280	4
10	0	0	210	4
11	0	0	140	4
12	0	0	70	4
13	0	0	0	4
14	0	0	910	4
15	0	0	980	4
16	3.9	0	1050	4
17	10	0	1085	4
18	20	0	1015	4
19	20	0	945	4
20	20	0	35	4
21	20	0	105	4
22	20	0	175	4
23	20	0	245	4
24	20	0	315	4
25	20	0	385	4
26	20	0	455	4
27	20	0	525	4
28	20	0	595	4
29	20	0	665	4
30	20	0	735	4
31	20	0	805	4
32	20	0	875	4
33	150	0	840	3
34	150	0	770	3
35	150	0	700	3
36	150	0	630	3
37	150	0	560	3
38	150	0	490	3
39	150	0	420	3
40	150	0	350	3
41	150	0	280	3
42	150	0	210	3
43	150	0	140	3
44	150	0	70	3
45	150	0	0	3
46	150	0	910	3
47	150	0	980	3
48	153.9	0	1050	3
49	160	0	1085	3
50	170	0	1015	3
51	170	0	945	3
52	170	0	35	3
53	170	0	105	3
54	170	0	175	3
55	170	0	245	3
56	170	0	315	3
57	170	0	385	3
58	170	0	455	3
59	170	0	525	3
60	170	0	595	3
61	170	0	665	3
62	170	0	735	3
63	170	0	805	3
64	170	0	875	3
65	300	0	840	2
66	300	0	770	2
67	300	0	700	2
68	300	0	630	2
69	300	0	560	2
70	300	0	490	2
71	300	0	420	2
72	300	0	350	2
73	300	0	280	2
74	300	0	210	2
75	300	0	140	2
76	300	0	70	2
77	300	0	0	2
78	300	0	910	2
79	300	0	980	2
80	303.9	0	1050	2
81	310	0	1085	2
82	320	0	1015	2
83	320	0	945	2
84	320	0	35	2
85	320	0	105	2
86	320	0	175	2
87	320	0	245	2
88	320	0	315	2
89	320	0	385	2
90	320	0	455	2
91	320	0	525	2
92	320	0	595	2
93	320	0	665	2
94	320	0	735	2
95	320	0	805	2
96	320	0	875	2
97	450	0	840	1
98	450	0	770	1
99	450	0	700	1
100	450	0	630	1
101	450	0	560	1
102	450	0	490	1
103	450	0	420	1
104	450	0	350	1
105	450	0	280	1
106	450	0	210	1
107	450	0	140	1
108	450	0	70	1
109	450	0	0	1
110	450	0	910	1
111	450	0	980	1
112	453.9	0	1050	1
113	460	0	1085	1
114	470	0	1015	1
115	470	0	945	1
116	470	0	35	1
117	470	0	105	1
118	470	0	175	1
119	470	0	245	1
120	470	0	315	1
121	470	0	385	1
122	470	0	455	1
123	470	0	525	1
124	470	0	595	1
125	470	0	665	1
126	470	0	735	1
127	470	0	805	1
128	470	0	875	1

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
for i=1:size(probewiring,1)
text(s.x(i)-5,s.z(i),num2str(s.channels(i)),'FontSize',9)
end
axis([min(s.x)-50 max(s.x)+50 min(s.z)-50 max(s.z)+50])
axis equal
set(gca,'FontSize',10,'TickDir','out')
