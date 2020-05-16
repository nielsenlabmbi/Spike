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

1	0	0	6200	2
2	0	0	6000	2
3	0	0	5800	2
4	0	0	5600	2
5	0	0	5400	2
6	0	0	5200	2
7	0	0	5000	2
8	0	0	4800	2
9	0	0	4600	2
10	0	0	4400	2
11	0	0	4200	2
12	0	0	4000	2
13	0	0	3800	2
14	0	0	3600	2
15	0	0	3400	2
16	0	0	3200	2
17	0	0	3000	2
18	0	0	2800	2
19	0	0	2600	2
20	0	0	2400	2
21	0	0	2200	2
22	0	0	2000	2
23	0	0	1800	2
24	0	0	1600	2
25	0	0	1400	2
26	0	0	1200	2
27	0	0	1000	2
28	0	0	800	2
29	0	0	600	2
30	0	0	400	2
31	0	0	200	2
32	0	0	0	2
33	3000	0	0	1
34	3000	0	200	1
35	3000	0	400	1
36	3000	0	600	1
37	3000	0	800	1
38	3000	0	1000	1
39	3000	0	1200	1
40	3000	0	1400	1
41	3000	0	1600	1
42	3000	0	1800	1
43	3000	0	2000	1
44	3000	0	2200	1
45	3000	0	2400	1
46	3000	0	2600	1
47	3000	0	2800	1
48	3000	0	3000	1
49	3000	0	3200	1
50	3000	0	3400	1
51	3000	0	3600	1
52	3000	0	3800	1
53	3000	0	4000	1
54	3000	0	4200	1
55	3000	0	4400	1
56	3000	0	4600	1
57	3000	0	4800	1
58	3000	0	5000	1
59	3000	0	5200	1
60	3000	0	5400	1
61	3000	0	5600	1
62	3000	0	5800	1
63	3000	0	6000	1
64	3000	0	6200	1

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
b=size(probewiring,2);disp(b);
for i=1:size(probewiring,1)
text(s.x(i)-5,s.z(i),num2str(s.channels(i)),'FontSize',9)
end
axis([min(s.x)-50 max(s.x)+50 min(s.z)-50 max(s.z)+50])
axis equal
set(gca,'FontSize',10,'TickDir','out')
