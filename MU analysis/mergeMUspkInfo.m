function mergeMUspkInfo(filepath,animalID,unitID,expID,probeID,startID,stopID)

%this function merges all MU spike files into one file for further analysis
%will remove events flagged as duplicates to stay consistent with the
%spkSort file
%
%input:
%filepath: path to the animal folder (e.g., Z:\EphysNew\Data)
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
%probeID: probe id
%startID: start job ID
%stopID: stop job ID
%
%output:
%MUspkMerge structure with fields
%spktimes: spike times
%detCh: detection channel for each spike
%detChSort: detection channels reorganized according to z position
%(bottom of probe is channel 1 in this organization)
%NDuplicate: number of cooccuring spike along the entire probe
%jobId: job file a spike is stored in

wb=waitbar(0,'Merging MUspkinfo');

expname=[animalID '_u' unitID '_' expID];

iniFile=0; %just in case the first file is empty
for i=startID:stopID
    waitbar((i-startID)/(stopID-startID),wb);
   
    %load file
    fname=fullfile(filepath,animalID,expname,'SpikeFiles',[expname '_j' num2str(i) '_p' num2str(probeID) '_MUspkinfo']);           
    load(fname); %generates spk
    
    %output
    if isfield(spk,'flagDuplicate') %necessary to guard against empty files
        if iniFile==0
            idx=find(spk.flagDuplicate==0);
            MUspkMerge.spktimes=spk.spkTimesDet(idx);
            MUspkMerge.detCh=spk.detCh(idx);
            MUspkMerge.detChSort=spk.detChSort(idx);
            MUspkMerge.NDuplicate=spk.NDuplicate (idx);
            MUspkMerge.jobId=ones(size(idx))*i;

            iniFile=1;
        else
            idx=find(spk.flagDuplicate==0);
            MUspkMerge.spktimes=[MUspkMerge.spktimes spk.spkTimesDet(idx)];
            MUspkMerge.detCh=[MUspkMerge.detCh spk.detCh(idx)];
            MUspkMerge.detChSort=[MUspkMerge.detChSort spk.detChSort(idx)];
            MUspkMerge.NDuplicate=[MUspkMerge.NDuplicate spk.NDuplicate(idx)];
            MUspkMerge.jobId=[MUspkMerge.jobId ones(size(idx))*i];
        end
        clear spk;
    end
end

save(fullfile(filepath,animalID,expname,[expname '_p' num2str(probeID) '_MUspkMerge']),'MUspkMerge'); 
close(wb);