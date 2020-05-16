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

1	0	0	600	2
2	0	0	550	2
3	0	0	500	2
4	0	0	450	2
5	0	0	400	2
6	0	0	350	2
7	0	0	300	2
8	0	0	250	2
9	0	0	200	2
10	0	0	150	2
11	0	0	100	2
12	0	0	50	2
13	0	0	0	2
14	0	0	650	2
15	0	0	700	2
16	3.9	0	750	2
17	10	0	775	2
18	20	0	725	2
19	20	0	675	2
20	20	0	25	2
21	20	0	75	2
22	20	0	125	2
23	20	0	175	2
24	20	0	225	2
25	20	0	275	2
26	20	0	325	2
27	20	0	375	2
28	20	0	425	2
29	20	0	475	2
30	20	0	525	2
31	20	0	575	2
32	20	0	625	2
33	300	0	600	1
34	300	0	550	1
35	300	0	500	1
36	300	0	450	1
37	300	0	400	1
38	300	0	350	1
39	300	0	300	1
40	300	0	250	1
41	300	0	200	1
42	300	0	150	1
43	300	0	100	1
44	300	0	50	1
45	300	0	0	1
46	300	0	650	1
47	300	0	700	1
48	304	0	750	1
49	310.1	0	775	1
50	320.1	0	725	1
51	320.1	0	675	1
52	320.1	0	25	1
53	320.1	0	75	1
54	320.1	0	125	1
55	320.1	0	175	1
56	320.1	0	225	1
57	320.1	0	275	1
58	320.1	0	325	1
59	320.1	0	375	1
60	320.1	0	425	1
61	320.1	0	475	1
62	320.1	0	525	1
63	320.1	0	575	1
64	320.1	0	625	1

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
