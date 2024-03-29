function resortSpkSort(expFolder,animalID,unitID,expID,probeID,copyToZ)
%to be able to use SortGui with the old format spike files, need to change
%the order spikes occur in spkSort (old format: sorted by time; new format:
%sorted by jobs first, then by channel, then time)
%do this after running the extractSpikesOldFormat function (needs the jobs
%boundaries)
% input parameters:
% expFolder - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% copyToZ - copy id file to Z



%load spkSort file
expname=[animalID '_u' unitID '_' expID];
load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort'])); %generates spkSort

%also load id file to get the job boundaries
load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); %generates id


%loop through jobs and get order
sortSpikes=[];
for j=1:length(id.extractSpikes.jobStart)
    %find the spikes in the job
    idx=find(spkSort.spktimes>=id.extractSpikes.jobStart(j) & spkSort.spktimes<id.extractSpikes.jobStop(j));
    
    %sort - channel first, spktimes after
    tmpArray=[spkSort.spktimes(idx);spkSort.detCh(idx)];
    [~,idSort]=sortrows(tmpArray',2);
    sortSpikes=[sortSpikes idx(idSort)];
end

%sort spikes now
spkSort.spktimes=spkSort.spktimes(sortSpikes);
spkSort.unitid=spkSort.unitid(sortSpikes);
spkSort.detCh=spkSort.detCh(sortSpikes);


%we also need to add some other info to spkSort for the sortGui to work
spkSort.info.jobStart=0;
spkSort.info.jobStop=length(id.extractSpikes.jobStart)-1;
spkSort.info.artRej='Off';

%save
save(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort']),'spkSort');

if copyToZ==1
    zbase='Z:\EphysNew\processedSpikes';
    save(fullfile(zbase,animalID,expname,[expname '_p' num2str(probeID) '_spkSort']),'spkSort');
end