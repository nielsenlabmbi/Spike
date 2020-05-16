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

1	0	0	50	2
2	0	0	0	2
3	0	0	99.9	2
4	0	0	150	2
5	0	0	200	2
6	0	0	250	2
7	0	0	300	2
8	0	0	350	2
9	0	0	400	2
10	4	0	450	2
11	20	0	475	2
12	20	0	375	2
13	20	0	275	2
14	20	0	175	2
15	20	0	75	2
16	20	0	-25	2
17	20	0	25	2
18	20	0	125	2
19	20	0	225	2
20	20	0	325	2
21	20	0	425	2
22	20	0	500	2
23	36	0	450	2
24	40	0	400	2
25	40	0	350	2
26	40	0	300	2
27	40	0	250	2
28	40	0	200	2
29	40	0	150	2
30	40	0	0	2
31	40	0	50	2
32	40	0	100	2
33	300.1	0	50	1
34	300.1	0	0	1
35	300.1	0	99.9	1
36	300.1	0	150	1
37	300.1	0	200	1
38	300.1	0	250	1
39	300.1	0	300	1
40	300.1	0	350	1
41	300.1	0	400	1
42	304.1	0	450	1
43	320.1	0	475	1
44	320.1	0	375	1
45	320.1	0	275	1
46	320.1	0	175	1
47	320.1	0	75	1
48	320.1	0	-25	1
49	320.1	0	25	1
50	320.1	0	125	1
51	320.1	0	225	1
52	320.1	0	325	1
53	320.1	0	425	1
54	320.1	0	500	1
55	336.1	0	450	1
56	340.1	0	400	1
57	340.1	0	350	1
58	340.1	0	300	1
59	340.1	0	250	1
60	340.1	0	200	1
61	340.1	0	150	1
62	340.1	0	0	1
63	340.1	0	50	1
64	340.1	0	100	1

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
