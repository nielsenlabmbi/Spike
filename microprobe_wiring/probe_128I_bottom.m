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
1	0	0	3000	4
2	0	0	2800	4
3	0	0	2600	4
4	0	0	2400	4
5	0	0	2200	4
6	0	0	2000	4
7	0	0	1800	4
8	0	0	1600	4
9	0	0	1400	4
10	0	0	1200	4
11	0	0	1000	4
12	0	0	800	4
13	0	0	600	4
14	0	0	400	4
15	0	0	200	4
16	0	0	0	4
17	0	0	100	4
18	0	0	300	4
19	0	0	500	4
20	0	0	700	4
21	0	0	900	4
22	0	0	1100	4
23	0	0	1300	4
24	0	0	1500	4
25	0	0	1700	4
26	0	0	1900	4
27	0	0	2100	4
28	0	0	2300	4
29	0	0	2500	4
30	0	0	2700	4
31	0	0	2900	4
32	0	0	3100	4
33	330	0	3000	3
34	330	0	2800	3
35	330	0	2600	3
36	330	0	2400	3
37	330	0	2200	3
38	330	0	2000	3
39	330	0	1800	3
40	330	0	1600	3
41	330	0	1400	3
42	330	0	1200	3
43	330	0	1000	3
44	330	0	800	3
45	330	0	600	3
46	330	0	400	3
47	330	0	200	3
48	330	0	0	3
49	330	0	100	3
50	330	0	300	3
51	330	0	500	3
52	330	0	700	3
53	330	0	900	3
54	330	0	1100	3
55	330	0	1300	3
56	330	0	1500	3
57	330	0	1700	3
58	330	0	1900	3
59	330	0	2100	3
60	330	0	2300	3
61	330	0	2500	3
62	330	0	2700	3
63	330	0	2900	3
64	330	0	3100	3
65	660	0	3000	2
66	660	0	2800	2
67	660	0	2600	2
68	660	0	2400	2
69	660	0	2200	2
70	660	0	2000	2
71	660	0	1800	2
72	660	0	1600	2
73	660	0	1400	2
74	660	0	1200	2
75	660	0	1000	2
76	660	0	800	2
77	660	0	600	2
78	660	0	400	2
79	660	0	200	2
80	660	0	0	2
81	660	0	100	2
82	660	0	300	2
83	660	0	500	2
84	660	0	700	2
85	660	0	900	2
86	660	0	1100	2
87	660	0	1300	2
88	660	0	1500	2
89	660	0	1700	2
90	660	0	1900	2
91	660	0	2100	2
92	660	0	2300	2
93	660	0	2500	2
94	660	0	2700	2
95	660	0	2900	2
96	660	0	3100	2
97	990	0	3000	1
98	990	0	2800	1
99	990	0	2600	1
100	990	0	2400	1
101	990	0	2200	1
102	990	0	2000	1
103	990	0	1800	1
104	990	0	1600	1
105	990	0	1400	1
106	990	0	1200	1
107	990	0	1000	1
108	990	0	800	1
109	990	0	600	1
110	990	0	400	1
111	990	0	200	1
112	990	0	0	1
113	990	0	100	1
114	990	0	300	1
115	990	0	500	1
116	990	0	700	1
117	990	0	900	1
118	990	0	1100	1
119	990	0	1300	1
120	990	0	1500	1
121	990	0	1700	1
122	990	0	1900	1
123	990	0	2100	1
124	990	0	2300	1
125	990	0	2500	1
126	990	0	2700	1
127	990	0	2900	1
128	990	0	3100	1

];


tipelectrode=25;    %nearest tip-electrode vertical distance in microns.

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
