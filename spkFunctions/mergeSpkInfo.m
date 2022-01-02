function mergeSpkInfo(filepath,expname,probeID,startID,stopID)

%this function merges all spike files into one file for MUA analysis
%it saves most of the same info that is contained in spkSort
%
%input:
%filepath: path to the animal folder (e.g.., Z:\EphysNew\Data\FEXX1)
%expname: name of experiment (FEXX1_u000_000)
%probeID: probe id
%startID: start job ID
%stopID: stop job ID

wb=waitbar(0,'Merging spkInfo');

for i=startID:stopID
    waitbar((i-startID)/(stopID-startID),wb);
   
    %load file
    fname=fullfile(filepath,expname,'SpikeFiles',[expname '_j' num2str(i) '_p' num2str(probeID) '_spkinfo']);           
    load(fname); %generates spk
    
    %output
    if i==startID
        spkMerge.spktimes=spk.spkTimesDet;
        spkMerge.detCh=spk.detCh;
        spkMerge.detChSort=spk.detChSort;
        spkMerge.jobId=ones(size(spk.spkTimesDet))*i;
    else
        spkMerge.spktimes=[spkMerge.spktimes spk.spkTimesDet];
        spkMerge.detCh=[spkMerge.detCh spk.detCh];
        spkMerge.detChSort=[spkMerge.detChSort spk.detChSort];
        spkMerge.jobId=[spkMerge.jobId ones(size(spk.spkTimesDet))*i];
    end
    clear spk;
end

save(fullfile(filepath,expname,[expname '_p' num2str(probeID) '_spkMrge']),'spkMerge'); 
close(wb);