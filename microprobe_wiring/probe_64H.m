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
1	0	0	190	2
2	0	0	163.75	2
3	0	0	137.5	2
4	0	0	110	2
5	0	0	82.5	2
6	0	0	55	2
7	0	0	27.5	2
8	0	0	0	2
9	0	0	217.5	2
10	2.5	0	245	2
11	22.5	0	286.25	2
12	22.5	0	258.75	2
13	22.5	0	231.25	2
14	22.5	0	203.75	2
15	22.5	0	-13.75	2
16	22.5	0	13.75	2
17	22.5	0	41.25	2
18	22.5	0	68.75	2
19	22.5	0	96.25	2
20	22.5	0	123.75	2
21	22.5	0	151.25	2
22	22.5	0	176.25	2
23	42.5	0	245	2
24	45	0	217.5	2
25	45	0	0	2
26	45	0	27.5	2
27	45	0	55	2
28	45	0	82.5	2
29	45	0	110	2
30	45	0	137.5	2
31	45	0	163.75	2
32	45	0	190	2
33	200.1	0	190	1
34	200.1	0	163.75	1
35	200.1	0	137.5	1
36	200.1	0	110	1
37	200.1	0	82.5	1
38	200.1	0	55	1
39	200.1	0	27.5	1
40	200.1	0	0	1
41	200.1	0	217.5	1
42	202.6	0	245	1
43	222.6	0	286.25	1
44	222.6	0	258.75	1
45	222.6	0	231.25	1
46	222.6	0	203.75	1
47	222.6	0	-13.75	1
48	222.6	0	13.75	1
49	222.6	0	41.25	1
50	222.6	0	68.75	1
51	222.6	0	96.25	1
52	222.6	0	123.75	1
53	222.6	0	151.25	1
54	222.6	0	176.25	1
55	242.6	0	245	1
56	245.1	0	217.5	1
57	245.1	0	0	1
58	245.1	0	27.5	1
59	245.1	0	55	1
60	245.1	0	82.5	1
61	245.1	0	110	1
62	245.1	0	137.5	1
63	245.1	0	163.75	1
64	245.1	0	190	1

];


tipelectrode=21;    %nearest tip-electrode vertical distance in microns.

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
