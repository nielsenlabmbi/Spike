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
1	0	0	700	2
2	0	0	650	2
3	0	0	600	2
4	0	0	550	2
5	0	0	500	2
6	0	0	450	2
7	0	0	400	2
8	0	0	350	2
9	0	0	300	2
10	0	0	250	2
11	0	0	200	2
12	0	0	150	2
13	0	0	100	2
14	0	0	50	2
15	0	0	0	2
16	0	0	750	2
17	0	0	800	2
18	0	0	850	2
19	0	0	900	2
20	0	0	950	2
21	4	0	1000	2
22	20	0	1050	2
23	20	0	975	2
24	20	0	875	2
25	20	0	775	2
26	20	0	675	2
27	20	0	575	2
28	20	0	475	2
29	20	0	375	2
30	20	0	275	2
31	20	0	175	2
32	20	0	75	2
33	20	0	25	2
34	20	0	125	2
35	20	0	225	2
36	20	0	325	2
37	20	0	425	2
38	20	0	525	2
39	20	0	625	2
40	20	0	725	2
41	20	0	825	2
42	20	0	925	2
43	20	0	1025	2
44	36	0	1000	2
45	40	0	950	2
46	40	0	900	2
47	40	0	850	2
48	40	0	800	2
49	40	0	0	2
50	40	0	50	2
51	40	0	100	2
52	40	0	150	2
53	40	0	200	2
54	40	0	250	2
55	40	0	300	2
56	40	0	350	2
57	40	0	400	2
58	40	0	450	2
59	40	0	500	2
60	40	0	550	2
61	40	0	600	2
62	40	0	650	2
63	40	0	700	2
64	40	0	750	2
65	400	0	700	1
66	400	0	650	1
67	400	0	600	1
68	400	0	550	1
69	400	0	500	1
70	400	0	450	1
71	400	0	400	1
72	400	0	350	1
73	400	0	300	1
74	400	0	250	1
75	400	0	200	1
76	400	0	150	1
77	400	0	100	1
78	400	0	50	1
79	400	0	0	1
80	400	0	750	1
81	400	0	800	1
82	400	0	850	1
83	400	0	900	1
84	400	0	950	1
85	404	0	1000	1
86	420	0	1050	1
87	420	0	975	1
88	420	0	875	1
89	420	0	775	1
90	420	0	675	1
91	420	0	575	1
92	420	0	475	1
93	420	0	375	1
94	420	0	275	1
95	420	0	175	1
96	420	0	75	1
97	420	0	25	1
98	420	0	125	1
99	420	0	225	1
100	420	0	325	1
101	420	0	425	1
102	420	0	525	1
103	420	0	625	1
104	420	0	725	1
105	420	0	825	1
106	420	0	925	1
107	420	0	1025	1
108	436	0	1000	1
109	440	0	950	1
110	440	0	900	1
111	440	0	850	1
112	440	0	800	1
113	440	0	0	1
114	440	0	50	1
115	440	0	100	1
116	440	0	150	1
117	440	0	200	1
118	440	0	250	1
119	440	0	300	1
120	440	0	350	1
121	440	0	400	1
122	440	0	450	1
123	440	0	500	1
124	440	0	550	1
125	440	0	600	1
126	440	0	650	1
127	440	0	700	1
128	440	0	750	1

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
