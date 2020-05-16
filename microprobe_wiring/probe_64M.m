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

1	0	0	1150	1
2	0	0	1100	1
3	0	0	1050	1
4	0	0	1000	1
5	0	0	950	1
6	0	0	900	1
7	0	0	850	1
8	0	0	800	1
9	0	0	750	1
10	0	0	700	1
11	0	0	650	1
12	0	0	600	1
13	0	0	550	1
14	0	0	500	1
15	0	0	450	1
16	0	0	400	1
17	0	0	350	1
18	0	0	300	1
19	0	0	250	1
20	0	0	200	1
21	0	0	150	1
22	0	0	100	1
23	0	0	50	1
24	0	0	0	1
25	0	0	1200	1
26	0	0	1250	1
27	0	0	1300	1
28	0	0	1350	1
29	0	0	1400	1
30	0	0	1450	1
31	0	0	1500	1
32	0	0	1550	1
33	25	0	1575	1
34	25	0	1525	1
35	25	0	1475	1
36	25	0	1425	1
37	25	0	1375	1
38	25	0	1325	1
39	25	0	1275	1
40	25	0	1225	1
41	25	0	1175	1
42	25	0	1125	1
43	25	0	1075	1
44	25	0	1025	1
45	25	0	975	1
46	25	0	925	1
47	25	0	875	1
48	25	0	825	1
49	25	0	775	1
50	25	0	725	1
51	25	0	675	1
52	25	0	625	1
53	25	0	575	1
54	25	0	525	1
55	25	0	475	1
56	25	0	425	1
57	25	0	375	1
58	25	0	325	1
59	25	0	275	1
60	25	0	225	1
61	25	0	175	1
62	25	0	125	1
63	25	0	75	1
64	25	0	25	1

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
