function convertProps(expFolder,animalID,unitID,expID,partialProp)

%transform new spike properties file into augustos format

% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% expFolder - base folder for experiments (string)
% partialProp - avoid for loop to speed things up


expname=[animalID '_u' unitID '_' expID];

load(fullfile(expFolder,animalID,expname,[expname  '_spk'])); %generates spkInfo

%output matrix
Properties=zeros(16,length(spkInfo.EnDet));

%minimal set of properties
Properties(5,:)=spkInfo.PksDet;
Properties(10,:)=spkInfo.AmpMinDet;
Properties(11,:)=spkInfo.EnDet;
Properties(12,:)=spkInfo.WidthDet;

Properties(14,:)=spkInfo.comZMinDet; %this is the closest we can get to Augustos center of mass value
Properties(15,:)=spkInfo.spkTimesDet;
Properties(16,:)=spkInfo.detCh;

%full set of properties
if partialProp~=1
   for i=1:length(spkInfo.EnAll)
       Amp=spkInfo.EnAll{i};
       %assumption: there will at least be 2 additional channels
       Properties(4,i)=Amp(2);
       Properties(6,i)=Amp(3);
       if length(Amp)>=5
           Properties(3,i)=Amp(4);
           Properties(7,i)=Amp(5);
       end
       if length(Amp)>=7
           Properties(2,i)=Amp(6);
           Properties(8,i)=Amp(7);
       end
       if length(Amp)>=9
           Properties(1,i)=Amp(8);
           Properties(9,i)=Amp(9);
       end
    
       Properties(13,i)=Amp(1)/Amp(end);
   end
       
end


idk = zeros(1,length(Properties(1,:)));
PropTitles{1} = 'Amp (-4)';
PropTitles{2} = 'Amp (-3)';
PropTitles{3} = 'Amp (-2)';
PropTitles{4} = 'Amp (-1)';
PropTitles{5} = 'Amp (0)';
PropTitles{6} = 'Amp (1)';
PropTitles{7} = 'Amp (2)';
PropTitles{8} = 'Amp (3)';
PropTitles{9} = 'Amp (4)';
PropTitles{10} = 'Pk2Pk Amp';
PropTitles{11} = 'Energy';
PropTitles{12} = 'Wvf width';
PropTitles{13} = 'CHs width';
PropTitles{14} = 'CHs pos';
PropTitles{15} = 'Time(samples)';
PropTitles{16} = 'Ch detected';


save(fullfile(expFolder,animalID,expname,[expname  '_Spikes.mat']),'Properties','idk','PropTitles')