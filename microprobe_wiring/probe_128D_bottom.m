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
33	330	0	600	3
34	330	0	550	3
35	330	0	500	3
36	330	0	450	3
37	330	0	400	3
38	330	0	350	3
39	330	0	300	3
40	330	0	250	3
41	330	0	200	3
42	330	0	150	3
43	330	0	100	3
44	330	0	50	3
45	330	0	0	3
46	330	0	650	3
47	330	0	700	3
48	333.9	0	750	3
49	340	0	775	3
50	350	0	725	3
51	350	0	675	3
52	350	0	25	3
53	350	0	75	3
54	350	0	125	3
55	350	0	175	3
56	350	0	225	3
57	350	0	275	3
58	350	0	325	3
59	350	0	375	3
60	350	0	425	3
61	350	0	475	3
62	350	0	525	3
63	350	0	575	3
64	350	0	625	3
65	660	0	600	2
66	660	0	550	2
67	660	0	500	2
68	660	0	450	2
69	660	0	400	2
70	660	0	350	2
71	660	0	300	2
72	660	0	250	2
73	660	0	200	2
74	660	0	150	2
75	660	0	100	2
76	660	0	50	2
77	660	0	0	2
78	660	0	650	2
79	660	0	700	2
80	663.9	0	750	2
81	670	0	775	2
82	680	0	725	2
83	680	0	675	2
84	680	0	25	2
85	680	0	75	2
86	680	0	125	2
87	680	0	175	2
88	680	0	225	2
89	680	0	275	2
90	680	0	325	2
91	680	0	375	2
92	680	0	425	2
93	680	0	475	2
94	680	0	525	2
95	680	0	575	2
96	680	0	625	2
97	990	0	600	1
98	990	0	550	1
99	990	0	500	1
100	990	0	450	1
101	990	0	400	1
102	990	0	350	1
103	990	0	300	1
104	990	0	250	1
105	990	0	200	1
106	990	0	150	1
107	990	0	100	1
108	990	0	50	1
109	990	0	0	1
110	990	0	650	1
111	990	0	700	1
112	993.9	0	750	1
113	1000	0	775	1
114	1010	0	725	1
115	1010	0	675	1
116	1010	0	25	1
117	1010	0	75	1
118	1010	0	125	1
119	1010	0	175	1
120	1010	0	225	1
121	1010	0	275	1
122	1010	0	325	1
123	1010	0	375	1
124	1010	0	425	1
125	1010	0	475	1
126	1010	0	525	1
127	1010	0	575	1
128	1010	0	625	1

];


tipelectrode=30;    %nearest tip-electrode vertical distance in microns.

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

% %To plot the labeled channels:
% figure(2)
% clf
% plot(s.x,s.z,'sqr', 'MarkerSize',11)
% hold on
% for i=1:size(probewiring,1)
% text(s.x(i)-5,s.z(i),num2str(s.channels(i)),'FontSize',9)
% end
% axis([min(s.x)-50 max(s.x)+50 min(s.z)-50 max(s.z)+50])
% axis equal
% set(gca,'FontSize',10,'TickDir','out')
