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

1	0	0	700	1
2	0	0	650	1
3	0	0	600	1
4	0	0	550	1
5	0	0	500	1
6	0	0	450	1
7	0	0	400	1
8	0	0	350	1
9	0	0	300	1
10	0	0	250	1
11	0	0	200	1
12	0	0	150	1
13	0	0	100	1
14	0	0	50	1
15	0	0	0	1
16	0	0	750	1
17	0	0	800	1
18	0	0	850	1
19	0	0	900	1
20	0	0	950	1
21	4	0	1000	1
22	20	0	1050	1
23	20	0	975	1
24	20	0	875	1
25	20	0	775	1
26	20	0	675	1
27	20	0	575	1
28	20	0	475	1
29	20	0	375	1
30	20	0	275	1
31	20	0	175	1
32	20	0	75	1
33	20	0	25	1
34	20	0	125	1
35	20	0	225	1
36	20	0	325	1
37	20	0	425	1
38	20	0	525	1
39	20	0	625	1
40	20	0	725	1
41	20	0	825	1
42	20	0	925	1
43	20	0	1025	1
44	36	0	1000	1
45	40	0	950	1
46	40	0	900	1
47	40	0	850	1
48	40	0	800	1
49	40	0	0	1
50	40	0	50	1
51	40	0	100	1
52	40	0	150	1
53	40	0	200	1
54	40	0	250	1
55	40	0	300	1
56	40	0	350	1
57	40	0	400	1
58	40	0	450	1
59	40	0	500	1
60	40	0	550	1
61	40	0	600	1
62	40	0	650	1
63	40	0	700	1
64	40	0	750	1

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
