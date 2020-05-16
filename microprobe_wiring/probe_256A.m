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

   1          40           0         750           4
           2          40           0         700           4
           3          40           0         650           4
           4          40           0         600           4
           5          40           0         550           4
           6          40           0         500           4
           7          40           0         450           4
           8          40           0         400           4
           9          40           0         350           4
          10          40           0         300           4
          11          40           0         250           4
          12          40           0         200           4
          13          40           0         150           4
          14          40           0         100           4
          15          40           0          50           4
          16          40           0           0           4
          17          40           0         800           4
          18          40           0         850           4
          19          40           0         900           4
          20          40           0         950           4
          21          36           0        1000           4
          22          20           0        1025           4
          23          20           0         925           4
          24          20           0         825           4
          25          20           0         725           4
          26          20           0         625           4
          27          20           0         525           4
          28          20           0         425           4
          29          20           0         325           4
          30          20           0         225           4
          31          20           0         125           4
          32          20           0          25           4
          33          20           0          75           4
          34          20           0         175           4
          35          20           0         275           4
          36          20           0         375           4
          37          20           0         475           4
          38          20           0         575           4
          39          20           0         675           4
          40          20           0         775           4
          41          20           0         875           4
          42          20           0         975           4
          43          20           0        1050           4
          44           4           0        1000           4
          45           0           0         950           4
          46           0           0         900           4
          47           0           0         850           4
          48           0           0         800           4
          49           0           0         750           4
          50           0           0           0           4
          51           0           0          50           4
          52           0           0         100           4
          53           0           0         150           4
          54           0           0         200           4
          55           0           0         250           4
          56           0           0         300           4
          57           0           0         350           4
          58           0           0         400           4
          59           0           0         450           4
          60           0           0         500           4
          61           0           0         550           4
          62           0           0         600           4
          63           0           0         650           4
          64           0           0         700           4
          65         400           0         700           3
          66         400           0         650           3
          67         400           0         600           3
          68         400           0         550           3
          69         400           0         500           3
          70         400           0         450           3
          71         400           0         400           3
          72         400           0         350           3
          73         400           0         300           3
          74         400           0         250           3
          75         400           0         200           3
          76         400           0         150           3
          77         400           0         100           3
          78         400           0          50           3
          79         400           0           0           3
          80         400           0         750           3
          81         400           0         800           3
          82         400           0         850           3
          83         400           0         900           3
          84         400           0         950           3
          85         404           0        1000           3
          86         420           0        1050           3
          87         420           0         975           3
          88         420           0         875           3
          89         420           0         775           3
          90         420           0         675           3
          91         420           0         575           3
          92         420           0         475           3
          93         420           0         375           3
          94         420           0         275           3
          95         420           0         175           3
          96         420           0          75           3
          97         420           0          25           3
          98         420           0         125           3
          99         420           0         225           3
         100         420           0         325           3
         101         420           0         425           3
         102         420           0         525           3
         103         420           0         625           3
         104         420           0         725           3
         105         420           0         825           3
         106         420           0         925           3
         107         420           0        1025           3
         108         436           0        1000           3
         109         440           0         950           3
         110         440           0         900           3
         111         440           0         850           3
         112         440           0         800           3
         113         440           0           0           3
         114         440           0          50           3
         115         440           0         100           3
         116         440           0         150           3
         117         440           0         200           3
         118         440           0         250           3
         119         440           0         300           3
         120         440           0         350           3
         121         440           0         400           3
         122         440           0         450           3
         123         440           0         500           3
         124         440           0         550           3
         125         440           0         600           3
         126         440           0         650           3
         127         440           0         700           3
         128         440           0         750           3
         129         800           0         700           2
         130         800           0         650           2
         131         800           0         600           2
         132         800           0         550           2
         133         800           0         500           2
         134         800           0         450           2
         135         800           0         400           2
         136         800           0         350           2
         137         800           0         300           2
         138         800           0         250           2
         139         800           0         200           2
         140         800           0         150           2
         141         800           0         100           2
         142         800           0          50           2
         143         800           0           0           2
         144         800           0         750           2
         145         800           0         800           2
         146         800           0         850           2
         147         800           0         900           2
         148         800           0         950           2
         149         804           0        1000           2
         150         820           0        1050           2
         151         820           0         975           2
         152         820           0         875           2
         153         820           0         775           2
         154         820           0         675           2
         155         820           0         575           2
         156         820           0         475           2
         157         820           0         375           2
         158         820           0         275           2
         159         820           0         175           2
         160         820           0          75           2
         161         820           0          25           2
         162         820           0         125           2
         163         820           0         225           2
         164         820           0         325           2
         165         820           0         425           2
         166         820           0         525           2
         167         820           0         625           2
         168         820           0         725           2
         169         820           0         825           2
         170         820           0         925           2
         171         820           0        1025           2
         172         836           0        1000           2
         173         840           0         950           2
         174         840           0         900           2
         175         840           0         850           2
         176         840           0         800           2
         177         840           0           0           2
         178         840           0          50           2
         179         840           0         100           2
         180         840           0         150           2
         181         840           0         200           2
         182         840           0         250           2
         183         840           0         300           2
         184         840           0         350           2
         185         840           0         400           2
         186         840           0         450           2
         187         840           0         500           2
         188         840           0         550           2
         189         840           0         600           2
         190         840           0         650           2
         191         840           0         700           2
         192         840           0         750           2
         193        1240           0         750           1
         194        1240           0         700           1
         195        1240           0         650           1
         196        1240           0         600           1
         197        1240           0         550           1
         198        1240           0         500           1
         199        1240           0         450           1
         200        1240           0         400           1
         201        1240           0         350           1
         202        1240           0         300           1
         203        1240           0         250           1
         204        1240           0         200           1
         205        1240           0         150           1
         206        1240           0         100           1
         207        1240           0          50           1
         208        1240           0           0           1
         209        1240           0         800           1
         210        1240           0         850           1
         211        1240           0         900           1
         212        1240           0         950           1
         213        1236           0        1000           1
         214        1220           0        1025           1
         215        1220           0         925           1
         216        1220           0         825           1
         217        1220           0         725           1
         218        1220           0         625           1
         219        1220           0         525           1
         220        1220           0         425           1
         221        1220           0         325           1
         222        1220           0         225           1
         223        1220           0         125           1
         224        1220           0          25           1
         225        1220           0          75           1
         226        1220           0         175           1
         227        1220           0         275           1
         228        1220           0         375           1
         229        1220           0         475           1
         230        1220           0         575           1
         231        1220           0         675           1
         232        1220           0         775           1
         233        1220           0         875           1
         234        1220           0         975           1
         235        1220           0        1050           1
         236        1204           0        1000           1
         237        1200           0         950           1
         238        1200           0         900           1
         239        1200           0         850           1
         240        1200           0         800           1
         241        1200           0         750           1
         242        1200           0           0           1
         243        1200           0          50           1
         244        1200           0         100           1
         245        1200           0         150           1
         246        1200           0         200           1
         247        1200           0         250           1
         248        1200           0         300           1
         249        1200           0         350           1
         250        1200           0         400           1
         251        1200           0         450           1
         252        1200           0         500           1
         253        1200           0         550           1
         254        1200           0         600           1
         255        1200           0         650           1
         256        1200           0         700           1
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
