% load([fileparts(mfilename('fullpath')) '/Settings.mat'])
% experiment = inputdlg('Enter experiment name: ');
% experiment = experiment{1};
load('Settings.mat')
f = uigetdir(expFolder,'Please select the experiment folder');
experiment = f(find(f == '/',1,'last')+1 : end);
expFolder = f(1:find(f == '/',1,'last')-1);            

PropertiesAll=[];
for i = 1:200
currentFile = [expFolder '/' experiment '/SpikeFiles/Spikes_' num2str(i-1) '.mat'];

load(currentFile)
PropertiesAll=[PropertiesAll Properties];
delete(currentFile)
end
Properties = PropertiesAll;
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

save([expFolder '/' experiment '/SpikeFiles/' experiment '_Spikes.mat'],'Properties','idk','PropTitles')