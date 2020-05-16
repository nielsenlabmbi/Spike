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
1	0	0	210	4
2	0	0	180	4
3	0	0	150	4
4	0	0	120	4
5	0	0	90	4
6	0	0	60	4
7	0	0	30	4
8	0	0	0	4
9	0	0	240	4
10	0	0	270	4
11	25	0	315	4
12	25	0	285	4
13	25	0	255	4
14	25	0	225	4
15	25	0	-15	4
16	25	0	15	4
17	25	0	45	4
18	25	0	75	4
19	25	0	105	4
20	25	0	135	4
21	25	0	165	4
22	25	0	195	4
23	50	0	270	4
24	50	0	240	4
25	50	0	0	4
26	50	0	30	4
27	50	0	60	4
28	50	0	90	4
29	50	0	120	4
30	50	0	150	4
31	50	0	180	4
32	50	0	210	4
33	200	0	210	3
34	200	0	180	3
35	200	0	150	3
36	200	0	120	3
37	200	0	90	3
38	200	0	60	3
39	200	0	30	3
40	200	0	0	3
41	200	0	240	3
42	200	0	270	3
43	225	0	315	3
44	225	0	285	3
45	225	0	255	3
46	225	0	225	3
47	225	0	-15	3
48	225	0	15	3
49	225	0	45	3
50	225	0	75	3
51	225	0	105	3
52	225	0	135	3
53	225	0	165	3
54	225	0	195	3
55	250	0	270	3
56	250	0	240	3
57	250	0	0	3
58	250	0	30	3
59	250	0	60	3
60	250	0	90	3
61	250	0	120	3
62	250	0	150	3
63	250	0	180	3
64	250	0	210	3
65	400	0	210	2
66	400	0	180	2
67	400	0	150	2
68	400	0	120	2
69	400	0	90	2
70	400	0	60	2
71	400	0	30	2
72	400	0	0	2
73	400	0	240	2
74	400	0	270	2
75	425	0	315	2
76	425	0	285	2
77	425	0	255	2
78	425	0	225	2
79	425	0	-15	2
80	425	0	15	2
81	425	0	45	2
82	425	0	75	2
83	425	0	105	2
84	425	0	135	2
85	425	0	165	2
86	425	0	195	2
87	450	0	270	2
88	450	0	240	2
89	450	0	0	2
90	450	0	30	2
91	450	0	60	2
92	450	0	90	2
93	450	0	120	2
94	450	0	150	2
95	450	0	180	2
96	450	0	210	2
97	600	0	210	1
98	600	0	180	1
99	600	0	150	1
100	600	0	120	1
101	600	0	90	1
102	600	0	60	1
103	600	0	30	1
104	600	0	0	1
105	600	0	240	1
106	600	0	270	1
107	625	0	315	1
108	625	0	285	1
109	625	0	255	1
110	625	0	225	1
111	625	0	-15	1
112	625	0	15	1
113	625	0	45	1
114	625	0	75	1
115	625	0	105	1
116	625	0	135	1
117	625	0	165	1
118	625	0	195	1
119	650	0	270	1
120	650	0	240	1
121	650	0	0	1
122	650	0	30	1
123	650	0	60	1
124	650	0	90	1
125	650	0	120	1
126	650	0	150	1
127	650	0	180	1
128	650	0	210	1

];


tipelectrode=35.5;    %nearest tip-electrode vertical distance in microns.

connector_position='bottom';

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
