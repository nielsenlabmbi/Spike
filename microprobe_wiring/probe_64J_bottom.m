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

1	0	0	6200	1
2	0	0	6000	1
3	0	0	5800	1
4	0	0	5600	1
5	0	0	5400	1
6	0	0	5200	1
7	0	0	5000	1
8	0	0	4800	1
9	0	0	4600	1
10	0	0	4400	1
11	0	0	4200	1
12	0	0	4000	1
13	0	0	3800	1
14	0	0	3600	1
15	0	0	3400	1
16	0	0	3200	1
17	0	0	3000	1
18	0	0	2800	1
19	0	0	2600	1
20	0	0	2400	1
21	0	0	2200	1
22	0	0	2000	1
23	0	0	1800	1
24	0	0	1600	1
25	0	0	1400	1
26	0	0	1200	1
27	0	0	1000	1
28	0	0	800	1
29	0	0	600	1
30	0	0	400	1
31	0	0	200	1
32	0	0	0	1
33	0	0	100	1
34	0	0	300	1
35	0	0	500	1
36	0	0	700	1
37	0	0	900	1
38	0	0	1100	1
39	0	0	1300	1
40	0	0	1500	1
41	0	0	1700	1
42	0	0	1900	1
43	0	0	2100	1
44	0	0	2300	1
45	0	0	2500	1
46	0	0	2700	1
47	0	0	2900	1
48	0	0	3100	1
49	0	0	3300	1
50	0	0	3500	1
51	0	0	3700	1
52	0	0	3900	1
53	0	0	4100	1
54	0	0	4300	1
55	0	0	4500	1
56	0	0	4700	1
57	0	0	4900	1
58	0	0	5100	1
59	0	0	5300	1
60	0	0	5500	1
61	0	0	5700	1
62	0	0	5900	1
63	0	0	6100	1
64	0	0	6300	1


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
