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
1	0	0	190	4
2	0	0	163.75	4
3	0	0	137.5	4
4	0	0	110	4
5	0	0	82.5	4
6	0	0	55	4
7	0	0	27.5	4
8	0	0	0	4
9	0	0	217.5	4
10	2.5	0	245	4
11	22.5	0	268.25	4
12	22.5	0	258.75	4
13	22.5	0	231.25	4
14	22.5	0	203.75	4
15	22.5	0	-13.75	4
16	22.5	0	13.75	4
17	22.5	0	41.25	4
18	22.5	0	68.75	4
19	22.5	0	96.25	4
20	22.5	0	123.75	4
21	22.5	0	151.25	4
22	22.5	0	176.25	4
23	42.5	0	245	4
24	45	0	217.5	4
25	45	0	0	4
26	45	0	27.5	4
27	45	0	55	4
28	45	0	82.5	4
29	45	0	110	4
30	45	0	137.5	4
31	45	0	163.75	4
32	45	0	190	4
33	150	0	190	3
34	150	0	163.75	3
35	150	0	137.5	3
36	150	0	110	3
37	150	0	82.5	3
38	150	0	55	3
39	150	0	27.5	3
40	150	0	0	3
41	150	0	217.5	3
42	152.25	0	245	3
43	172.5	0	268.25	3
44	172.5	0	258.75	3
45	172.5	0	231.25	3
46	172.5	0	203.75	3
47	172.5	0	-13.75	3
48	172.5	0	13.75	3
49	172.5	0	41.25	3
50	172.5	0	68.75	3
51	172.5	0	96.25	3
52	172.5	0	123.75	3
53	172.5	0	151.25	3
54	172.5	0	176.25	3
55	192.5	0	245	3
56	195	0	217.5	3
57	195	0	0	3
58	195	0	27.5	3
59	195	0	55	3
60	195	0	82.5	3
61	195	0	110	3
62	195	0	137.5	3
63	195	0	163.75	3
64	195	0	190	3
65	300	0	190	2
66	300	0	163.75	2
67	300	0	137.5	2
68	300	0	110	2
69	300	0	82.5	2
70	300	0	55	2
71	300	0	27.5	2
72	300	0	0	2
73	300	0	217.5	2
74	302.5	0	245	2
75	322.5	0	268.25	2
76	322.5	0	258.75	2
77	322.5	0	231.25	2
78	322.5	0	203.75	2
79	322.5	0	-13.75	2
80	322.5	0	13.75	2
81	322.5	0	41.25	2
82	322.5	0	68.75	2
83	322.5	0	96.25	2
84	322.5	0	123.75	2
85	322.5	0	151.25	2
86	322.5	0	176.25	2
87	342.5	0	245	2
88	345	0	217.5	2
89	345	0	0	2
90	345	0	27.5	2
91	345	0	55	2
92	345	0	82.5	2
93	345	0	110	2
94	345	0	137.5	2
95	345	0	163.75	2
96	345	0	190	2
97	450	0	190	1
98	450	0	163.75	1
99	450	0	137.5	1
100	450	0	110	1
101	450	0	82.5	1
102	450	0	55	1
103	450	0	27.5	1
104	450	0	0	1
105	450	0	217.5	1
106	452.25	0	245	1
107	472.5	0	268.25	1
108	472.5	0	258.75	1
109	472.5	0	231.25	1
110	472.5	0	203.75	1
111	472.5	0	-13.75	1
112	472.5	0	13.75	1
113	472.5	0	41.25	1
114	472.5	0	68.75	1
115	472.5	0	96.25	1
116	472.5	0	123.75	1
117	472.5	0	151.25	1
118	472.5	0	176.25	1
119	492.5	0	245	1
120	495	0	217.5	1
121	495	0	0	1
122	495	0	27.5	1
123	495	0	55	1
124	495	0	82.5	1
125	495	0	110	1
126	495	0	137.5	1
127	495	0	163.75	1
128	495	0	190	1

];


tipelectrode=20.5;    %nearest tip-electrode vertical distance in microns.

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
